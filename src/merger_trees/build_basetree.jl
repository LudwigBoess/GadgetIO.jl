"""
Port of L-BaseTree (`MergerTrees/L-BaseTree/main.c`).

Given two consecutive subfind catalogues `A` (source) and `B` (target),
the descendant of each subhalo in `A` is the subhalo in `B` that contains
the largest *rank-weighted* fraction of `A`'s particles:

    score(B_j | A_i) = Σ_{p ∈ A_i ∩ B_j}  1 / (rank_A(p) + 1)^α

with `α = 2/3` and `rank_A(p)` being the (0-based) position of particle
`p` within the most-bound list of `A_i`. The descendant is the `B_j` with
the largest score. This rank weighting favours the most-bound (typically
inner) particles, which are robust under stripping.
"""

const BASETREE_ALPHA = 2.0 / 3.0

"""
    SubhaloCatalogue

Thin in-memory view of the per-subhalo data needed for descendant
matching: particle-ID list and per-subhalo offsets/lengths.

| Field        | Meaning                                                  |
| :----------- | :------------------------------------------------------- |
| `sublen`     | number of particle IDs per subhalo                       |
| `suboffset`  | start of each subhalo within `ids` (0-based)             |
| `ids`        | concatenated particle IDs in most-bound order            |
"""
struct SubhaloCatalogue
    sublen::Vector{Int32}
    suboffset::Vector{Int64}
    ids::Vector{UInt64}
end

Base.length(c::SubhaloCatalogue) = length(c.sublen)

"""
    load_subhalo_catalogue(sub_base) -> SubhaloCatalogue

Read SLEN, SOFF and PID from a subfind output. SLEN/SOFF and PID
typically live on *different* subfind parttypes (subhalo-type vs
ID-type in LGADGET3), so we issue one batched read per parttype.

For multi-file subfind outputs, different L-Gadget/SUBFIND flavours
write SOFF differently:

* `NEW_SUBFIND` (without `LGADGET3`) stores **global** SOFF already.
* `LGADGET3` (and the legacy non-`NEW_SUBFIND` path) store SOFF
  **local** to each sub-file's PID block and require an `idcount`
  bump — `cat->SubOffset[subcount + i] = suboffset_temp[i] + idcount`
  in `L-BaseTree/main.c`.

Rather than introduce a compile-time switch, this function picks the
representation that produces a consistent `(suboffset + sublen)` ≤
`length(ids)` layout. The check is cheap (one pass over sublen) and
robust as long as the catalogue isn't truly malformed.
"""
function load_subhalo_catalogue(sub_base::AbstractString)
    parttype_sub = subfind_block_parttype(sub_base, "SLEN")
    parttype_pid = subfind_block_parttype(sub_base, "PID")

    sub_blocks = read_subfind_blocks(sub_base, ["SLEN", "SOFF"]; parttype = parttype_sub)
    pid_blocks = read_subfind_blocks(sub_base, ["PID"];          parttype = parttype_pid)

    sublen        = Int32.(sub_blocks["SLEN"])
    suboffset_raw = Int64.(sub_blocks["SOFF"])
    ids           = UInt64.(pid_blocks["PID"])

    suboffset = _resolve_suboffsets(sub_base, suboffset_raw, sublen,
                                    Int64(length(ids)),
                                    parttype_sub, parttype_pid)

    return SubhaloCatalogue(sublen, suboffset, ids)
end

# Choose between "SOFF as written" (already-global flavour) and the
# inter-sub-file bumped version. Picks the first one that produces
# a layout consistent with the PID list length.
function _resolve_suboffsets(sub_base::AbstractString,
                             raw::Vector{Int64},
                             sublen::Vector{Int32},
                             n_ids::Int64,
                             parttype_sub::Integer,
                             parttype_pid::Integer)

    isempty(sublen) && return raw

    raw_max = _max_suboffset_end(raw, sublen)
    if raw_max <= n_ids
        # SOFF is already valid as a global index into the concatenated
        # PID list (NEW_SUBFIND flavour). Use as-is.
        return raw
    end

    # Otherwise treat SOFF as per-sub-file local and bump.
    corrected = _global_suboffsets(sub_base, raw, parttype_sub, parttype_pid)
    corr_max  = _max_suboffset_end(corrected, sublen)
    if corr_max <= n_ids
        return corrected
    end

    throw(ErrorException(
        "SubhaloCatalogue layout inconsistent: max (suboffset + sublen) is " *
        "$raw_max as written and $corr_max after the inter-sub-file " *
        "correction, but length(ids) = $n_ids. This usually means PID is " *
        "on an unexpected subfind parttype, or the SOFF/SLEN/PID blocks " *
        "come from incompatible sub-files."))
end

@inline function _max_suboffset_end(suboffset::Vector{Int64},
                                    sublen::Vector{Int32})
    m = Int64(0)
    @inbounds for i in eachindex(sublen)
        e = suboffset[i] + Int64(sublen[i])
        m = max(m, e)
    end
    return m
end

# Bump per-sub-file local SOFF values into the global IdList frame.
function _global_suboffsets(sub_base::AbstractString,
                            suboffset_local::Vector{Int64},
                            parttype_sub::Integer,
                            parttype_pid::Integer)
    h_outer = read_header(sub_base)
    num_files = isfile(sub_base) ? 1 : Int(h_outer.num_files)
    num_files == 1 && return suboffset_local

    suboffset = copy(suboffset_local)
    sub_start = 0          # 0-based subhalo index of first subhalo in current file
    pid_start = Int64(0)   # cumulative PID count from earlier files
    @inbounds for fi in 1:num_files
        h = read_header(select_file(sub_base, fi - 1))
        n_sub_fi = Int(h.npart[parttype_sub + 1])
        n_pid_fi = Int64(h.npart[parttype_pid + 1])
        for i in 1:n_sub_fi
            suboffset[sub_start + i] += pid_start
        end
        sub_start += n_sub_fi
        pid_start += n_pid_fi
    end
    return suboffset
end


"""
    IdToHaloMap

Sorted parallel-array lookup mapping a particle ID to its (0-based)
subhalo index in some target catalogue. Replaces the
`Dict{UInt64,Int32}` used in the first port — ~3× less memory, ~2-3×
faster lookups under load, and trivially shareable across threads
(read-only after construction).
"""
struct IdToHaloMap
    ids::Vector{UInt64}    # sorted ascending
    halos::Vector{Int32}   # parallel
end

Base.length(m::IdToHaloMap) = length(m.ids)

"""
    build_id_to_halo(cat::SubhaloCatalogue) -> IdToHaloMap

Map every particle ID owned by some subhalo back to that subhalo's
0-based index. Sorted by ID for `searchsortedfirst` lookup.
"""
function build_id_to_halo(cat::SubhaloCatalogue)
    n = isempty(cat.sublen) ? 0 : Int(sum(cat.sublen))

    # one-time bounds sanity check — without this, a SOFF that wasn't
    # globalised across sub-files segfaults inside the `@inbounds` read
    # below instead of producing a clear error.
    if n > 0
        max_end = Int64(0)
        @inbounds for i in eachindex(cat.sublen)
            e = cat.suboffset[i] + Int64(cat.sublen[i])
            max_end = max(max_end, e)
        end
        max_end > length(cat.ids) && throw(ErrorException(
            "SubhaloCatalogue layout inconsistent: max (suboffset + sublen) = " *
            "$max_end exceeds length(ids) = $(length(cat.ids)). For multi-file " *
            "subfind outputs, SOFF must be globalised across sub-files."))
    end

    ids_flat   = Vector{UInt64}(undef, n)
    halos_flat = Vector{Int32}(undef, n)

    k = 0
    @inbounds for i in eachindex(cat.sublen)
        off = cat.suboffset[i]
        len = Int(cat.sublen[i])
        h = Int32(i - 1)   # 0-based subhalo index
        for j in 1:len
            k += 1
            ids_flat[k]   = cat.ids[off + j]
            halos_flat[k] = h
        end
    end

    perm = sortperm(ids_flat)
    return IdToHaloMap(ids_flat[perm], halos_flat[perm])
end

@inline function lookup(m::IdToHaloMap, id::UInt64)
    n = length(m.ids)
    n == 0 && return NO_HALO
    @inbounds begin
        idx = searchsortedfirst(m.ids, id)
        (idx > n || m.ids[idx] != id) && return NO_HALO
        return m.halos[idx]
    end
end


"""
    determine_descendants(catA, catB; alpha=2/3, snapnum_B=-1)
        -> (descendant_haloindex, descendant_snapnum, descendant_weight)

For each subhalo in `catA`, find the subhalo in `catB` that maximises
the rank-weighted ID-overlap score. Returns three vectors of length
`length(catA)`:

* `descendant_haloindex[i]` — 0-based descendant index in `catB`, `-1` if none
* `descendant_snapnum[i]`   — `snapnum_B` if found, `-1` otherwise
* `descendant_weight[i]`    — score of the chosen descendant

Mirrors the C `determine_descendants` in L-BaseTree: collect per-halo
candidate (haloB, weight) pairs, sort by `haloB`, sweep to aggregate
weights per haloB, and pick the argmax. Uses per-task scratch buffers
under `Threads.@spawn`-based parallelism so it is safe across Julia
versions including the post-1.10 task-migration model.
"""
function determine_descendants(catA::SubhaloCatalogue, catB::SubhaloCatalogue;
                               alpha::Real = BASETREE_ALPHA,
                               snapnum_B::Integer = -1)

    id_to_halo_B = build_id_to_halo(catB)

    nA = length(catA)
    descendant_haloindex = fill(NO_HALO, nA)
    descendant_snapnum   = fill(Int32(-1), nA)
    descendant_weight    = zeros(Float32, nA)

    nA == 0 && return descendant_haloindex, descendant_snapnum, descendant_weight

    maxlen = Int(maximum(catA.sublen))
    snap_B = Int32(snapnum_B)

    # parallel chunks; one Task per chunk, each with its own scratch
    nt = max(1, Threads.nthreads())
    chunks = Vector{Int}(undef, nt + 1)
    @inbounds for t in 0:nt
        chunks[t + 1] = round(Int, t * nA / nt)
    end

    @sync for t in 1:nt
        i_lo, i_hi = chunks[t] + 1, chunks[t + 1]
        i_lo > i_hi && continue
        Threads.@spawn _determine_descendants_chunk!(
            descendant_haloindex, descendant_snapnum, descendant_weight,
            catA, id_to_halo_B, maxlen, Float32(alpha), snap_B,
            i_lo, i_hi)
    end

    return descendant_haloindex, descendant_snapnum, descendant_weight
end

function _determine_descendants_chunk!(desc_halo::Vector{Int32},
                                       desc_snap::Vector{Int32},
                                       desc_w::Vector{Float32},
                                       catA::SubhaloCatalogue,
                                       id_to_halo_B::IdToHaloMap,
                                       maxlen::Int,
                                       alpha::Float32,
                                       snap_B::Int32,
                                       i_lo::Int, i_hi::Int)

    # per-task scratch — never shared, immune to task migration
    candlist = Vector{Tuple{Int32,Float32}}(undef, maxlen)

    @inbounds for i in i_lo:i_hi
        off = catA.suboffset[i]
        len = Int(catA.sublen[i])
        ncand = 0

        for j in 1:len
            id = catA.ids[off + j]
            haloB = lookup(id_to_halo_B, id)
            haloB == NO_HALO && continue
            ncand += 1
            w = Float32(1.0 / (Float32(j) ^ alpha))
            candlist[ncand] = (haloB, w)
        end

        ncand == 0 && continue

        # sort by haloindex, then sweep + pick max in one pass
        sort!(view(candlist, 1:ncand); by = first)

        prev_h  = Int32(-2)
        acc_w   = 0f0
        best_h  = NO_HALO
        best_w  = 0f0
        for j in 1:ncand
            h, w = candlist[j]
            if h == prev_h
                acc_w += w
            else
                if acc_w > best_w
                    best_w = acc_w
                    best_h = prev_h
                end
                prev_h = h
                acc_w  = w
            end
        end
        # flush last run
        if acc_w > best_w
            best_w = acc_w
            best_h = prev_h
        end

        desc_halo[i] = best_h
        desc_snap[i] = snap_B
        desc_w[i]    = best_w
    end

    return nothing
end


"""
    count_progenitors(descendant_haloindex, n_target) -> Vector{Int32}

Count, for each subhalo in the target snapshot, how many source subhalos
point to it as their descendant.
"""
function count_progenitors(descendant_haloindex::AbstractVector{<:Integer},
                           n_target::Integer)
    counts = zeros(Int32, n_target)
    for d in descendant_haloindex
        d >= 0 && (counts[d + 1] += Int32(1))
    end
    return counts
end


"""
    decide_upon_descendant!(descA_B, descA_C, descSnapA, weightA_B, weightA_C,
                            countB, countC; weight_fak=3.0)

Re-route subhalos in snap A from their primary descendant in B to the
secondary one in C when one of the following holds:

1. The B-descendant has multiple progenitors **and** the C-descendant
   has none — i.e. the halo seems to merge in B but survive in C.
2. The C-route weight (skipping snap B) is more than `weight_fak`× the
   B-route weight (only if `skip_by_weight=true`).
3. No B-descendant exists but a C-descendant does.

This is L-BaseTree's `decide_upon_descendant` ported with explicit
arguments (no globals). Modifies the descendant arrays in place.
"""
function decide_upon_descendant!(descA_B::AbstractVector{Int32},
                                 descA_C::AbstractVector{Int32},
                                 descSnapA_B::AbstractVector{Int32},
                                 descSnapA_C::AbstractVector{Int32},
                                 weightA_B::AbstractVector{Float32},
                                 weightA_C::AbstractVector{Float32},
                                 countB::AbstractVector{Int32},
                                 countC::AbstractVector{Int32};
                                 weight_fak::Real = 3.0,
                                 skip_by_weight::Bool = false)

    @assert length(descA_B) == length(descA_C) == length(descSnapA_B) ==
            length(descSnapA_C) == length(weightA_B) == length(weightA_C)

    n_rerouted_merge   = 0
    n_rerouted_weight  = 0
    n_rerouted_nob     = 0

    @inbounds for i in eachindex(descA_B)
        ib = descA_B[i]
        ic = descA_C[i]

        if ib >= 0 && ic >= 0
            # case 1: descendant in B is part of a merger, descendant in C is fresh
            if countB[ib + 1] > 1 && countC[ic + 1] == 0
                countB[ib + 1]   -= Int32(1)
                countC[ic + 1]   += Int32(1)
                descA_B[i]        = ic
                descSnapA_B[i]    = descSnapA_C[i]
                n_rerouted_merge += 1

            # case 2: optional weight-based override
            elseif skip_by_weight && weightA_C[i] / weight_fak > weightA_B[i]
                countB[ib + 1]    -= Int32(1)
                countC[ic + 1]    += Int32(1)
                descA_B[i]         = ic
                descSnapA_B[i]     = descSnapA_C[i]
                n_rerouted_weight += 1
            end

        elseif ib < 0 && ic >= 0
            # case 3: no immediate descendant but one survives skipping a snap
            descA_B[i]      = ic
            descSnapA_B[i]  = descSnapA_C[i]
            countC[ic + 1] += Int32(1)
            n_rerouted_nob += 1
        end
    end

    return (; n_rerouted_merge, n_rerouted_weight, n_rerouted_nob)
end


"""
    build_basetree(sub_base_A, sub_base_B, snapnum_B;
                   sub_base_C=nothing, snapnum_C=nothing,
                   alpha=2/3, skip_by_weight=false, weight_fak=3.0)
        -> SubDesc

End-to-end port of one L-BaseTree run for source snap A:

1. Load subhalo catalogues for A and B (and optionally C).
2. Find descendants A → B (primary).
3. If C is provided, find descendants A → C and apply the
   `decide_upon_descendant` rerouting.

Returns a `SubDesc` ready to be written with `write_sub_desc`.
"""
function build_basetree(sub_base_A::AbstractString,
                        sub_base_B::AbstractString,
                        snapnum_B::Integer;
                        sub_base_C::Union{Nothing,AbstractString} = nothing,
                        snapnum_C::Union{Nothing,Integer}         = nothing,
                        alpha::Real = BASETREE_ALPHA,
                        skip_by_weight::Bool = false,
                        weight_fak::Real = 3.0,
                        verbose::Bool = true)

    verbose && @info "Loading catalogues" sub_base_A sub_base_B sub_base_C
    catA = load_subhalo_catalogue(sub_base_A)
    catB = load_subhalo_catalogue(sub_base_B)

    descA_B, snapA_B, weightA_B =
        determine_descendants(catA, catB; alpha, snapnum_B)

    if sub_base_C !== nothing && snapnum_C !== nothing
        catC = load_subhalo_catalogue(sub_base_C)

        # secondary A→C and intermediate B→C
        descA_C, snapA_C, weightA_C =
            determine_descendants(catA, catC; alpha, snapnum_B = snapnum_C)
        descB_C, _, _ =
            determine_descendants(catB, catC; alpha, snapnum_B = snapnum_C)

        countB = count_progenitors(descA_B, length(catB))
        countC = count_progenitors(descB_C, length(catC))

        decide_upon_descendant!(descA_B, descA_C, snapA_B, snapA_C,
                                 weightA_B, weightA_C, countB, countC;
                                 weight_fak, skip_by_weight)
    end

    return SubDesc(Int32(length(catA)), descA_B, snapA_B, Int64[])
end
