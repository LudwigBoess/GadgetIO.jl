"""
Port of L-HaloTrees (`MergerTrees/L-HaloTrees/main.c`).

Combines per-snapshot subfind catalogues and `sub_desc_NNN` descendant
files (produced by L-BaseTree / `build_basetree`) into a global `Halo`
array, links FOF groups and progenitor chains, then walks each root in
depth-first order to emit a `MergerTreeFile`.

The pointer convention is 0-based (`-1` = none) on disk; conversion to
1-based Julia indices happens only inside the helpers exposed for users
(`walk_tree`, layout/visualisation).
"""

"""
    HaloAux

Scratch flags used while walking the trees. Mirrors L-HaloTrees'
`halo_aux_data` C struct.
"""
mutable struct HaloAux
    used::Bool
    fof_done::Bool
    target_index::Int32
    origin::Int32
end
HaloAux() = HaloAux(false, false, Int32(0), Int32(0))


"""
    LinkTable

Struct-of-Arrays scratch for the linked-list / DFS phase of tree
construction. Holds only the six fields that are touched during
`set_progenitor_pointers!` and `walk_it!` (5 link pointers + `Len`),
24 bytes per halo instead of the 104-byte `MergerTreeHalo`. This:

* lets `set_progenitor_pointers!` update a single `Int32` per link
  instead of copying a 104-byte immutable struct,
* keeps the hot-fields working set ~4× more cache-dense during DFS.

`FirstProgenitor` and `NextProgenitor` start as `NO_HALO` and are filled
in by `set_progenitor_pointers!`.
"""
struct LinkTable
    descendant::Vector{Int32}
    first_progenitor::Vector{Int32}
    next_progenitor::Vector{Int32}
    first_halo_in_fof::Vector{Int32}
    next_halo_in_fof::Vector{Int32}
    len::Vector{Int32}
end

LinkTable(n::Integer) = LinkTable(
    fill(NO_HALO, n), fill(NO_HALO, n), fill(NO_HALO, n),
    fill(NO_HALO, n), fill(NO_HALO, n), zeros(Int32, n))

Base.length(l::LinkTable) = length(l.descendant)


"""
    SubfindMeta

Container for the per-snapshot subfind data needed to build the global
halo array. Provide one entry per snapshot. All vectors are flat over
*all subfind files* of that snapshot, in subfind index order.
"""
struct SubfindMeta
    sublen::Vector{Int32}
    nsub_per_halo::Vector{Int32}    # length = number of FOF groups
    halo_M_Mean200::Vector{Float32}
    halo_M_Crit200::Vector{Float32}
    halo_M_TopHat::Vector{Float32}
    subpos::Matrix{Float32}         # 3 × nsub
    subvel::Matrix{Float32}         # 3 × nsub
    subspin::Matrix{Float32}        # 3 × nsub
    subveldisp::Vector{Float32}
    subvmax::Vector{Float32}
    subhalfmass::Vector{Float32}
    subMostBoundID::Vector{Int64}
    file_nr_of_sub::Vector{Int32}
    subhalo_index::Vector{Int32}
end

"""
    load_subfind_meta(sub_base; mass1="MVIR", mass2="M5CC", mass3="M500") -> SubfindMeta

Load the subfind blocks required for tree construction from a multi-file
or single-file subfind output. The three mass strings configure which
overdensity definitions go into `halo_M_TopHat`, `halo_M_Crit200`, and
`halo_M_Mean200` respectively (matching the C macros `MASS1_STRING`,
`MASS2_STRING`, `MASS3_STRING`).
"""
function load_subfind_meta(sub_base::AbstractString;
                           mass1::AbstractString = "MVIR",
                           mass2::AbstractString = "M5CC",
                           mass3::AbstractString = "M500")

    # split the requested blocks by their subfind parttype (FOF-halo vs
    # subhalo), then issue ONE batched read per parttype rather than 14
    # separate full-pipeline reads. This is the dominant runtime saving
    # on multi-file subfind outputs (Lustre/NFS).
    halo_blocks = [mass1, mass2, mass3, "NSUB"]
    sub_blocks  = ["SLEN", "SPOS", "SVEL", "SPIN",
                   "DSUB", "VMAX", "RHMS", "MBID"]

    halo_pt = subfind_block_parttype(sub_base, "NSUB")
    sub_pt  = subfind_block_parttype(sub_base, "SLEN")

    halo_data = read_subfind_blocks(sub_base, halo_blocks; parttype = halo_pt)
    sub_data  = read_subfind_blocks(sub_base, sub_blocks;  parttype = sub_pt)

    sublen          = Int32.(sub_data["SLEN"])
    nsub            = length(sublen)
    nsub_per_halo   = Int32.(halo_data["NSUB"])

    halo_M_TopHat  = Float32.(halo_data[mass1])
    halo_M_Crit200 = Float32.(halo_data[mass2])
    halo_M_Mean200 = Float32.(halo_data[mass3])

    subpos      = Float32.(sub_data["SPOS"])
    subvel      = Float32.(sub_data["SVEL"])
    subspin     = Float32.(sub_data["SPIN"])
    subveldisp  = Float32.(sub_data["DSUB"])
    subvmax     = Float32.(sub_data["VMAX"])
    subhalfmass = Float32.(sub_data["RHMS"])
    subMostBoundID = Int64.(sub_data["MBID"])

    # file_nr / subhalo_index are bookkeeping fields the C code records
    # while iterating subfind sub-files. Without that loop they can stay
    # at the trivial single-file values (filenr=0, idx=sub-index).
    file_nr_of_sub = zeros(Int32, nsub)
    subhalo_index  = Int32.(collect(0:nsub-1))

    return SubfindMeta(sublen, nsub_per_halo,
                       halo_M_Mean200, halo_M_Crit200, halo_M_TopHat,
                       subpos, subvel, subspin,
                       subveldisp, subvmax, subhalfmass,
                       subMostBoundID,
                       file_nr_of_sub, subhalo_index)
end


"""
    assemble_global_halos(metas, descendants, first_snap, last_snap)
        -> (halos::Vector{MergerTreeHalo}, links::LinkTable, first_in_snap::Dict{Int,Int32})

Build the global halo array from per-snapshot subfind meta and
descendant files. Halos are concatenated in reverse-snapshot order
(`last_snap`, `last_snap-1`, …) which matches the C code so that root
halos come first.

`metas[snap]` is a `SubfindMeta`. `descendants[snap]` is a `SubDesc`
(only required for `snap < last_snap`). The returned `links::LinkTable`
holds the SoA pointer + `Len` fields used by `set_progenitor_pointers!`
and `walk_it!`. `first_in_snap[snap]` returns the 0-based offset of the
snapshot's first halo in the global array.

The `MergerTreeHalo` entries are written with `NO_HALO` placeholders in
their `FirstProgenitor`/`NextProgenitor` slots; the final values are
merged in by `generate_trees` when it emits per-tree halos.
"""
function assemble_global_halos(metas::AbstractDict{Int,SubfindMeta},
                               descendants::AbstractDict{Int,SubDesc},
                               first_snap::Integer,
                               last_snap::Integer)

    tot_halos = sum(length(metas[s].sublen) for s in first_snap:last_snap)
    halos = Vector{MergerTreeHalo}(undef, tot_halos)
    links = LinkTable(tot_halos)

    first_in_snap = Dict{Int,Int32}()

    count = 0
    for num in last_snap:-1:first_snap
        first_in_snap[num] = Int32(count)
        meta = metas[num]

        if num < last_snap
            desc = descendants[num]
        else
            desc = nothing
        end

        nh = count
        sc = 0   # 0-based running subhalo index within this snapshot
        for gr in 1:length(meta.nsub_per_halo)
            gr_nh_start = nh
            nsub_in_grp = Int(meta.nsub_per_halo[gr])

            for subgr in 1:nsub_in_grp
                # FOF linkage
                first_halo_in_fof = Int32(gr_nh_start)
                next_halo_in_fof  = (subgr == nsub_in_grp) ? NO_HALO : Int32(nh + 1)

                # global descendant index
                if desc !== nothing
                    dh = desc.descendant_haloindex[sc + 1]
                    if dh >= 0
                        dsnap = Int(desc.descendant_snapnum[sc + 1])
                        descendant = first_in_snap[dsnap] + dh
                    else
                        descendant = NO_HALO
                    end
                else
                    descendant = NO_HALO
                end

                if subgr == 1
                    m200_mean = meta.halo_M_Mean200[gr]
                    m200_crit = meta.halo_M_Crit200[gr]
                    m_tophat  = meta.halo_M_TopHat[gr]
                else
                    m200_mean = 0f0; m200_crit = 0f0; m_tophat = 0f0
                end

                halos[nh + 1] = MergerTreeHalo(
                    descendant,
                    NO_HALO,   # FirstProgenitor — placeholder; merged in by generate_trees
                    NO_HALO,   # NextProgenitor  — placeholder; merged in by generate_trees
                    first_halo_in_fof,
                    next_halo_in_fof,
                    meta.sublen[sc + 1],
                    m200_mean, m200_crit, m_tophat,
                    (meta.subpos[1, sc + 1], meta.subpos[2, sc + 1], meta.subpos[3, sc + 1]),
                    (meta.subvel[1, sc + 1], meta.subvel[2, sc + 1], meta.subvel[3, sc + 1]),
                    meta.subveldisp[sc + 1],
                    meta.subvmax[sc + 1],
                    (meta.subspin[1, sc + 1], meta.subspin[2, sc + 1], meta.subspin[3, sc + 1]),
                    meta.subMostBoundID[sc + 1],
                    Int32(num),
                    meta.file_nr_of_sub[sc + 1],
                    meta.subhalo_index[sc + 1],
                    meta.subhalfmass[sc + 1],
                )

                # parallel SoA pointer fields, ~25× cheaper to update later
                links.descendant[nh + 1]        = descendant
                links.first_halo_in_fof[nh + 1] = first_halo_in_fof
                links.next_halo_in_fof[nh + 1]  = next_halo_in_fof
                links.len[nh + 1]               = meta.sublen[sc + 1]

                sc += 1
                nh += 1
            end
        end

        count = nh
    end

    return halos, links, first_in_snap
end


"""
    set_progenitor_pointers!(links::LinkTable)

Build the `FirstProgenitor`/`NextProgenitor` linked lists from the
already-filled `Descendant` field. Mirrors L-HaloTrees'
`set_progenitor_pointers`: at each descendant, progenitors are inserted
into a linked list so that the **largest** halo (by `Len`) sits at
`FirstProgenitor`. Within the chain, ordering matches the C code so
that downstream consumers (galaxybox-style traversals) see the same
sequence.

Operates on the SoA `LinkTable` so each insertion is a single `Int32`
write rather than a 104-byte struct copy.
"""
function set_progenitor_pointers!(links::LinkTable)
    n = length(links)
    @inbounds for i in 1:n
        desc = Int(links.descendant[i])
        desc < 0 && continue

        first = Int(links.first_progenitor[desc + 1])

        if first >= 0
            if links.len[i] >= links.len[first + 1]
                links.next_progenitor[i]        = Int32(first)
                links.first_progenitor[desc+1]  = Int32(i - 1)
            else
                links.next_progenitor[i]        = links.next_progenitor[first + 1]
                links.next_progenitor[first+1]  = Int32(i - 1)
            end
        else
            links.first_progenitor[desc + 1] = Int32(i - 1)
        end
    end
    return links
end

"""
    set_progenitor_pointers!(halos::Vector{MergerTreeHalo})

Convenience method that operates directly on the immutable halo array
by extracting the link fields into a `LinkTable`, linking, then
merging the new `FirstProgenitor` / `NextProgenitor` values back into
`halos`. Useful for callers that don't go through the full
`build_halotree` pipeline (e.g. tests that build a global array by
hand). For the pipeline itself prefer the `LinkTable` overload — it
saves the final materialisation pass.
"""
function set_progenitor_pointers!(halos::Vector{MergerTreeHalo})
    n = length(halos)
    links = LinkTable(n)
    @inbounds for i in 1:n
        h = halos[i]
        links.descendant[i]        = h.Descendant
        links.first_progenitor[i]  = h.FirstProgenitor
        links.next_progenitor[i]   = h.NextProgenitor
        links.first_halo_in_fof[i] = h.FirstHaloInFOFgroup
        links.next_halo_in_fof[i]  = h.NextHaloInFOFgroup
        links.len[i]               = h.Len
    end
    set_progenitor_pointers!(links)
    @inbounds for i in 1:n
        h = halos[i]
        halos[i] = MergerTreeHalo(
            h.Descendant,
            links.first_progenitor[i],
            links.next_progenitor[i],
            h.FirstHaloInFOFgroup, h.NextHaloInFOFgroup,
            h.Len, h.M_Mean200, h.M_Crit200, h.M_TopHat,
            h.Pos, h.Vel, h.VelDisp, h.Vmax, h.Spin,
            h.MostBoundID, h.SnapNum, h.FileNr, h.SubhaloIndex, h.SubhalfMass)
    end
    return halos
end


"""
    walk_it!(links, aux, root, callback)

Depth-first traversal used by L-HaloTrees to flatten a tree into a
contiguous depth-first array. `callback(i_global, k_local)` is invoked
for each visited halo (`i_global` is the 0-based global index, `k_local`
is its 0-based position inside the current tree). Operates on the
SoA `LinkTable` for cache density — only 24 bytes per halo touched
during the walk.
"""
function walk_it!(links::LinkTable,
                  aux::Vector{HaloAux},
                  root::Int32,
                  callback::F) where {F}
    counter = Ref(Int32(0))
    _walk_it!(links, aux, root, counter, callback)
    return counter[]
end

function _walk_it!(links::LinkTable,
                   aux::Vector{HaloAux},
                   i::Int32,
                   counter::Ref{Int32},
                   callback::F) where {F}
    @inbounds aux[i + 1].used = true
    @inbounds aux[i + 1].target_index = counter[]
    callback(i, counter[])
    counter[] += Int32(1)

    @inbounds desc = links.descendant[i + 1]
    if desc >= 0 && !(@inbounds aux[desc + 1].used)
        _walk_it!(links, aux, desc, counter, callback)
    end

    @inbounds p = links.first_progenitor[i + 1]
    while p >= 0
        if !(@inbounds aux[p + 1].used)
            _walk_it!(links, aux, p, counter, callback)
        end
        @inbounds p = links.next_progenitor[p + 1]
    end

    @inbounds fof = links.first_halo_in_fof[i + 1]
    if fof >= 0 && !(@inbounds aux[fof + 1].fof_done)
        @inbounds aux[fof + 1].fof_done = true
        q = fof
        while q >= 0
            if !(@inbounds aux[q + 1].used)
                _walk_it!(links, aux, q, counter, callback)
            end
            @inbounds q = links.next_halo_in_fof[q + 1]
        end
    end
end


"""
    generate_trees(halos, links, last_snap_nsub) -> MergerTreeFile

Walk each subhalo of the final snapshot that has not been visited yet
and flatten its connected tree depth-first. Pointer fields are
remapped from global halo indices to per-tree indices so the resulting
`MergerTreeFile` is self-contained.

`last_snap_nsub` is the number of subhalos in the final (root)
snapshot, i.e. `length(metas[last_snap].sublen)`. These are the
candidate roots — see L-HaloTrees `generate_trees`.

Single-pass: one DFS per tree records both `aux.target_index` and the
per-tree `origin` (global indices in DFS order). The per-tree
`MergerTreeHalo` output is constructed directly from `halos` + `links`
with pointers remapped via `aux.target_index` — no intermediate
verbatim-copy traversal.
"""
function generate_trees(halos::Vector{MergerTreeHalo},
                        links::LinkTable,
                        last_snap_nsub::Integer)

    n_total = length(halos)
    @assert length(links) == n_total
    aux = [HaloAux() for _ in 1:n_total]

    trees  = MergerTree[]
    origin = Int32[]   # reused per-tree scratch

    for i in 0:last_snap_nsub - 1
        @inbounds aux[i + 1].used && continue

        empty!(origin)
        walk_it!(links, aux, Int32(i),
                 (g, _k) -> push!(origin, g))

        n = length(origin)
        out = Vector{MergerTreeHalo}(undef, n)
        @inbounds for k in 1:n
            g = Int(origin[k])
            h = halos[g + 1]
            out[k] = MergerTreeHalo(
                _remap(aux, h.Descendant),
                _remap(aux, links.first_progenitor[g + 1]),
                _remap(aux, links.next_progenitor[g + 1]),
                _remap(aux, h.FirstHaloInFOFgroup),
                _remap(aux, h.NextHaloInFOFgroup),
                h.Len, h.M_Mean200, h.M_Crit200, h.M_TopHat,
                h.Pos, h.Vel, h.VelDisp, h.Vmax, h.Spin,
                h.MostBoundID, h.SnapNum, h.FileNr, h.SubhaloIndex, h.SubhalfMass,
            )
        end

        push!(trees, MergerTree(out))
    end

    tot_used = Int32(sum(length(t) for t in trees; init = 0))
    return MergerTreeFile(trees, Int32(length(trees)), tot_used)
end

@inline _remap(aux::Vector{HaloAux}, p::Int32) =
    p < 0 ? NO_HALO : (@inbounds aux[p + 1].target_index)


"""
    generate_trees(halos::Vector{MergerTreeHalo}, last_snap_nsub) -> MergerTreeFile

Convenience overload: extract a `LinkTable` from `halos` and dispatch
to the SoA implementation. Useful for tests and callers that build the
global halo array by hand with `FirstProgenitor` / `NextProgenitor`
already filled in.
"""
function generate_trees(halos::Vector{MergerTreeHalo},
                        last_snap_nsub::Integer)
    n = length(halos)
    links = LinkTable(n)
    @inbounds for i in 1:n
        h = halos[i]
        links.descendant[i]        = h.Descendant
        links.first_progenitor[i]  = h.FirstProgenitor
        links.next_progenitor[i]   = h.NextProgenitor
        links.first_halo_in_fof[i] = h.FirstHaloInFOFgroup
        links.next_halo_in_fof[i]  = h.NextHaloInFOFgroup
        links.len[i]               = h.Len
    end
    return generate_trees(halos, links, last_snap_nsub)
end


"""
    build_halotree(metas::Dict{Int,SubfindMeta},
                   descendants::Dict{Int,SubDesc},
                   first_snap, last_snap) -> MergerTreeFile

End-to-end port of L-HaloTrees: assemble the global halo array, set
progenitor pointers, then flatten into a `MergerTreeFile` ready to be
written with `write_merger_tree_file`.
"""
function build_halotree(metas::AbstractDict{Int,SubfindMeta},
                        descendants::AbstractDict{Int,SubDesc},
                        first_snap::Integer,
                        last_snap::Integer)

    halos, links, _ = assemble_global_halos(metas, descendants,
                                            first_snap, last_snap)
    set_progenitor_pointers!(links)
    return generate_trees(halos, links, length(metas[last_snap].sublen))
end
