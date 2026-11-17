"""
Lagrangian tracking of (sub)halo cores by their member particle IDs, plus a
particle-ID descendant linker.

Halo finders and merger-tree builders frequently lose an infalling halo during a
close passage (consistent-trees cleans it, subfind's unbinding strips it). The
particle IDs, however, are persistent across snapshots, so following a fixed set
of a halo's innermost ("most-bound-like") particles recovers a smooth trajectory
and a robust descendant link even where the trees fail.

Typical workflow:

```julia
# 1) seed: pick halo centres/radii at a snapshot where they are distinct
tags = tag_innermost_ids(snap_base_seed, centers, radii; parttype = 1, n_tag = 2000)

# 2) track every core across all snapshots (one read per snapshot)
snap_bases = Dict(s => "…/snapdir_\$(s)/snap_\$(s)" for s in 25:45)
tracks = track_id_groups(snap_bases, tags; parttype = 1)
#   tracks[g][snap] = (center, nfound, rms)

# 3) (optional) re-link two snapshots' halos by shared particle IDs
desc, nshared = link_by_shared_ids(tags_snapA, tags_snapB)
```

All positions are in the snapshot's native code units; these routines are
unit-agnostic.
"""

"""
    shrinking_sphere_center(pos; shrink=0.9, fmin=0.1, nmin=20) -> Vector{Float64}

Iterative shrinking-sphere centre of the 3×N position matrix `pos` (columns =
particles). Starts at the centre of mass and repeatedly keeps the particles
within `shrink`×(current maximum radius) of the running centre, until fewer than
`max(nmin, fmin*N)` particles remain (or the set stops shrinking). Robust to the
tagged cloud being partially stripped or captured.
"""
function shrinking_sphere_center(pos::AbstractMatrix{<:Real};
                                 shrink::Real = 0.9, fmin::Real = 0.1, nmin::Integer = 20)
    N = size(pos, 2)
    N == 0 && return fill(NaN, 3)
    idx = collect(1:N)
    cx = sum(@view pos[1, :]) / N
    cy = sum(@view pos[2, :]) / N
    cz = sum(@view pos[3, :]) / N
    target = max(nmin, ceil(Int, fmin * N))
    while length(idx) > target
        rmax2 = 0.0
        @inbounds for i in idx
            r2 = (pos[1, i] - cx)^2 + (pos[2, i] - cy)^2 + (pos[3, i] - cz)^2
            r2 > rmax2 && (rmax2 = r2)
        end
        thr2 = (shrink * shrink) * rmax2
        keep = Int[]
        @inbounds for i in idx
            ((pos[1, i] - cx)^2 + (pos[2, i] - cy)^2 + (pos[3, i] - cz)^2) <= thr2 && push!(keep, i)
        end
        (isempty(keep) || length(keep) == length(idx)) && break
        idx = keep
        n = length(idx)
        sx = sy = sz = 0.0
        @inbounds for i in idx
            sx += pos[1, i]; sy += pos[2, i]; sz += pos[3, i]
        end
        cx = sx / n; cy = sy / n; cz = sz / n
    end
    return [cx, cy, cz]
end


"""
    tag_innermost_ids(snap_base, centers, radii; parttype=1, n_tag=2000, verbose=false)
        -> Vector{Vector}

Read particle type `parttype` (`POS` and `ID`) from the (possibly distributed)
snapshot `snap_base` once and, for each halo `i` with centre `centers[i]`
(3-vector) and search radius `radii[i]`, return the IDs of up to `n_tag`
particles closest to the centre (within `radii[i]`). These "core tags" seed
[`track_id_groups`](@ref) / [`track_ids`](@ref) and are the linking unit for
[`link_by_shared_ids`](@ref). Positions and radii are in code units.
"""
function tag_innermost_ids(snap_base::AbstractString,
                           centers::AbstractVector, radii::AbstractVector;
                           parttype::Integer = 1, n_tag::Integer = 2000,
                           verbose::Bool = false)
    verbose && @info "tag_innermost_ids: reading POS/ID (type $parttype) from $snap_base"
    pos = read_block(snap_base, "POS"; parttype = parttype)
    ids = read_block(snap_base, "ID";  parttype = parttype)

    rad2  = Float64.(radii) .^ 2
    nc    = length(centers)
    bucket = [Tuple{Float64,Int}[] for _ in 1:nc]   # (r², particle index) within each radius

    @inbounds for i in axes(pos, 2)
        px = pos[1, i]; py = pos[2, i]; pz = pos[3, i]
        for h in 1:nc
            c = centers[h]
            r2 = (px - c[1])^2 + (py - c[2])^2 + (pz - c[3])^2
            r2 <= rad2[h] && push!(bucket[h], (r2, i))
        end
    end

    tags = Vector{Vector{eltype(ids)}}(undef, nc)
    for h in 1:nc
        b = bucket[h]
        sort!(b, by = first)
        nt = min(n_tag, length(b))
        tags[h] = eltype(ids)[ids[b[k][2]] for k in 1:nt]
        verbose && @info "  halo $h: tagged $nt of $(length(b)) within radius"
    end
    return tags
end


# centre / nfound / rms of a matched position cloud
@inline function _cloud_stats(pos, shrink, fmin, nmin)
    n = size(pos, 2)
    n == 0 && return (center = fill(NaN, 3), nfound = 0, rms = NaN)
    c = shrinking_sphere_center(pos; shrink, fmin, nmin)
    s = 0.0
    @inbounds for i in 1:n
        s += (pos[1, i] - c[1])^2 + (pos[2, i] - c[2])^2 + (pos[3, i] - c[3])^2
    end
    return (center = c, nfound = n, rms = sqrt(s / n))
end


"""
    track_ids(snap_bases, ids; parttype=1, verbose=false, shrink=0.9, fmin=0.1, nmin=20,
              pos0=nothing, r0=0.0, use_keys=false) -> Dict{Int,NamedTuple}

Follow a single Lagrangian set of particle `ids` (type `parttype`) across the
snapshots in `snap_bases` (a `snapnum => snap_base` map). For each snapshot reads
those IDs' positions and returns `(center, nfound, rms)`: the shrinking-sphere
centre, the number of IDs located, and the rms distance to the centre (a size /
stripping proxy). See [`track_id_groups`](@ref) to follow many sets at once.

`pos0`/`r0`/`use_keys` are forwarded to [`read_particles_by_id`](@ref): give a
search centre+radius (and `use_keys=true`) to use Peano-Hilbert keys for fast
region reads when a `.key` index is available.
"""
function track_ids(snap_bases::AbstractDict{<:Integer,<:AbstractString},
                   ids::AbstractVector{<:Integer};
                   parttype::Integer = 1, verbose::Bool = false,
                   shrink::Real = 0.9, fmin::Real = 0.1, nmin::Integer = 20,
                   pos0 = nothing, r0::Real = 0.0, use_keys::Bool = false)
    return track_id_groups(snap_bases, [collect(ids)];
                           parttype, verbose, shrink, fmin, nmin,
                           pos0, r0, use_keys)[1]
end


"""
    track_id_groups(snap_bases, groups; parttype=1, verbose=false,
                    shrink=0.9, fmin=0.1, nmin=20) -> Vector{Dict{Int,NamedTuple}}

Follow many Lagrangian sets at once. `groups` is a vector of particle-ID vectors
(e.g. the output of [`tag_innermost_ids`](@ref)). For each snapshot the **union**
of all IDs is read once with [`read_particles_by_id`](@ref), the tagged particles
are bucketed back to their groups, and every group's shrinking-sphere centre is
computed from that read — so the cost is one read per snapshot, not one per
group. Returns one `snapnum => (center, nfound, rms)` map per group. Groups are
assumed ~disjoint (a shared ID is attributed to the last group listing it).

`pos0`/`r0`/`use_keys` are forwarded to [`read_particles_by_id`](@ref); pass a
search centre+radius and `use_keys=true` to use Peano-Hilbert keys for fast
region reads when a `.key` index exists (otherwise the whole ID block is
scanned, which still works).
"""
function track_id_groups(snap_bases::AbstractDict{<:Integer,<:AbstractString},
                         groups::AbstractVector;
                         parttype::Integer = 1, verbose::Bool = false,
                         shrink::Real = 0.9, fmin::Real = 0.1, nmin::Integer = 20,
                         pos0 = nothing, r0::Real = 0.0, use_keys::Bool = false)
    union_ids = unique(reduce(vcat, groups))
    id2g = Dict{eltype(union_ids),Int}()
    for (g, grp) in enumerate(groups), x in grp
        id2g[x] = g
    end
    res = [Dict{Int,NamedTuple{(:center, :nfound, :rms),Tuple{Vector{Float64},Int,Float64}}}()
           for _ in groups]
    for snap in sort(collect(keys(snap_bases)))
        verbose && @info "track_id_groups: reading snap $snap"
        data = read_particles_by_id(snap_bases[snap], union_ids, ["POS"];
                                    parttype, verbose, pos0, r0, use_keys)
        pos = data["POS"]; ids = data["ID"]
        buckets = [Int[] for _ in groups]
        @inbounds for k in eachindex(ids)
            g = get(id2g, ids[k], 0)
            g > 0 && push!(buckets[g], k)
        end
        for g in eachindex(groups)
            sub = isempty(buckets[g]) ? zeros(eltype(pos), 3, 0) : pos[:, buckets[g]]
            res[g][snap] = _cloud_stats(sub, shrink, fmin, nmin)
        end
    end
    return res
end


"""
    merger_trees_from_descendants(keys, snaps, masses, desc_keys;
                                  no_descendant = -1, particle_mass = 1.0,
                                  positions = nothing, mostbound = nothing) -> Vector{MergerTree}

Build a forest of [`MergerTree`](@ref)s from a flat descendant table — e.g. the
particle-ID re-link produced by `link_by_shared_ids` (`relink_rockstar.jl`), or
any external catalogue. Each node `i` has a unique `keys[i]`, a `snaps[i]`
(snapshot number → `SnapNum`), a `masses[i]` (→ `M_TopHat`, and `Len =
masses/particle_mass` so the most massive progenitor heads each chain), and
`desc_keys[i]` = the key of its descendant (or `no_descendant` if it is a root).

Optionally pass `positions` (a vector of 3-tuples/vectors, → `Pos`) and
`mostbound` (a vector of ids, → `MostBoundID`) per node so the resulting trees
carry positions (e.g. for 3D trajectories) and a back-reference to the source
catalogue.

Nodes are grouped into connected trees by following the descendant chain to its
root; within each tree the pointers are made local (0-based) and
`FirstProgenitor`/`NextProgenitor` are filled by [`set_progenitor_pointers!`](@ref),
so the result plots directly with `MakieHelper.plot_merger_forest`.

```julia
trees = merger_trees_from_descendants(keys, snaps, mvir, desc_keys;
                                      no_descendant = (-1, -1), particle_mass = 3.87e7)
plot_merger_forest(trees; snap_to_scale = s -> a_of[s], mass_of = h -> log10(h.M_TopHat))
```
"""
function merger_trees_from_descendants(keys::AbstractVector, snaps::AbstractVector,
                                       masses::AbstractVector, desc_keys::AbstractVector;
                                       no_descendant = -1, particle_mass::Real = 1.0,
                                       positions = nothing, mostbound = nothing)
    n = length(keys)
    key2idx = Dict(keys[i] => i for i in 1:n)        # node key -> global 1-based index
    gdesc = fill(-1, n)                              # global 0-based descendant, -1 = root
    for i in 1:n
        dk = desc_keys[i]
        if dk != no_descendant && haskey(key2idx, dk) && key2idx[dk] != i
            gdesc[i] = key2idx[dk] - 1
        end
    end

    # assign each node to its root by following the descendant chain (memoised)
    rootof = zeros(Int, n)
    function findroot(i)
        rootof[i] != 0 && return rootof[i]
        d = gdesc[i]
        rootof[i] = d < 0 ? i : findroot(d + 1)
        return rootof[i]
    end
    for i in 1:n
        findroot(i)
    end

    groups = Dict{Int,Vector{Int}}()
    for i in 1:n
        push!(get!(groups, rootof[i], Int[]), i)
    end

    maxlen = Float64(typemax(Int32))
    trees = MergerTree[]
    for members in values(groups)
        g2l = Dict(members[k] => k for k in eachindex(members))    # global -> local 1-based
        halos = Vector{MergerTreeHalo}(undef, length(members))
        for k in eachindex(members)
            gi = members[k]
            ldesc = gdesc[gi] < 0 ? NO_HALO : Int32(g2l[gdesc[gi] + 1] - 1)
            len = Int32(clamp(round(Float64(masses[gi]) / particle_mass), 1.0, maxlen))
            pos = positions === nothing ? (0f0, 0f0, 0f0) :
                  (Float32(positions[gi][1]), Float32(positions[gi][2]), Float32(positions[gi][3]))
            mb  = mostbound === nothing ? Int64(0) : Int64(mostbound[gi])
            halos[k] = MergerTreeHalo(
                ldesc, NO_HALO, NO_HALO, NO_HALO, NO_HALO,
                len, 0f0, 0f0, Float32(masses[gi]),
                pos, (0f0, 0f0, 0f0), 0f0, 0f0, (0f0, 0f0, 0f0),
                mb, Int32(snaps[gi]), Int32(0), Int32(0), 0f0)
        end
        set_progenitor_pointers!(halos)
        push!(trees, MergerTree(halos))
    end
    return trees
end


"""
    link_by_shared_ids(tags_a, tags_b; min_shared=1) -> (descendant, nshared)

L-BaseTree-style descendant matching by shared particle IDs. `tags_a` and
`tags_b` are vectors of particle-ID sets for the halos of snapshot A and B
(e.g. from [`tag_innermost_ids`](@ref)). For each halo `i` in A, returns the
index into B of the halo sharing the most of A's IDs (`descendant[i]`, `0` if
fewer than `min_shared`) and that shared count (`nshared[i]`). Build a merger
tree by calling this on consecutive snapshots.
"""
function link_by_shared_ids(tags_a::AbstractVector, tags_b::AbstractVector;
                            min_shared::Integer = 1)
    setsb = [Set(t) for t in tags_b]
    descendant = zeros(Int, length(tags_a))
    nshared    = zeros(Int, length(tags_a))
    for (i, ta) in enumerate(tags_a)
        best = 0; bestn = 0
        for (j, sb) in enumerate(setsb)
            cnt = 0
            for x in ta
                (x in sb) && (cnt += 1)
            end
            cnt > bestn && (bestn = cnt; best = j)
        end
        if bestn >= min_shared
            descendant[i] = best
            nshared[i]    = bestn
        end
    end
    return descendant, nshared
end
