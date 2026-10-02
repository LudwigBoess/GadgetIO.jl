"""
    counting_sort(perm, digit, n)

Stable reorder of the particle indices `perm` by `digit[perm]`, with digits in `0:n-1`.
"""
function counting_sort(perm, digit, n)

    start = zeros(Int, n + 1)
    @inbounds for i ∈ perm
        start[digit[i]+2] += 1
    end
    for k = 2:n+1
        start[k] += start[k-1]
    end

    sorted = similar(perm)
    @inbounds for i ∈ perm
        k = digit[i] + 1
        start[k] += 1
        sorted[start[k]] = i
    end

    return sorted
end

"""
    find_root!(parent, i)

Root of `i` in the union-find forest `parent`, halving the path on the way.
"""
@inline function find_root!(parent, i)
    @inbounds while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end

"""
    link_cells!(parent, tree_size, a, b, x, cell_first, l2)

Join the groups of cells `a` and `b` if any particle pair between them is closer than
`√l2`. Particles of cell `a` are `x[:, cell_first[a]:cell_first[a+1]-1]`.
"""
@inline function link_cells!(parent, tree_size, a, b, x, cell_first, l2)

    ra = find_root!(parent, a)
    rb = find_root!(parent, b)
    ra == rb && return

    @inbounds for i = cell_first[a]:cell_first[a+1]-1
        xi, yi, zi = Float64(x[1, i]), Float64(x[2, i]), Float64(x[3, i])
        for j = cell_first[b]:cell_first[b+1]-1
            if (xi - x[1, j])^2 + (yi - x[2, j])^2 + (zi - x[3, j])^2 <= l2
                # union by size
                if tree_size[ra] < tree_size[rb]
                    ra, rb = rb, ra
                end
                parent[rb] = ra
                tree_size[ra] += tree_size[rb]
                return
            end
        end
    end
end

"""
    fof_labels(pos::AbstractMatrix{<:Real}, linking_length::Real)

Friends-of-friends group label of every particle in `pos` (3×N). Particles share a label if
a chain of pairs no further apart than `linking_length` connects them.

Particles are sorted into cells of side `0.999 linking_length/√3`, so all particles of a cell
are linked and the union-find runs over cells. Only cells up to two cells apart along each
axis can hold a linked pair. Every such cell pair is visited once and skipped if both cells
are in the same group already, which is the common case in dense regions.

The cells are indexed by a dense `nx × ny` array of columns, which is only small for a
compact particle set such as the high-resolution region of a zoom simulation.
"""
function fof_labels(pos::AbstractMatrix{<:Real}, linking_length::Real)

    N = size(pos, 2)
    iszero(N) && return Int32[]
    linking_length > 0 || error("The linking length must be positive, got $linking_length.")

    l2 = Float64(linking_length)^2
    c = 0.999 * linking_length / √3

    lo = Float64.(vec(minimum(pos, dims=2)))
    extent = (Float64.(vec(maximum(pos, dims=2))) .- lo) ./ c
    if maximum(extent) > 2^24 || extent[1] * extent[2] > 2^28
        error(@sprintf("Particles span %.3g × %.3g × %.3g cells of the linking length. ", extent...) *
              "fof_labels is meant for compact particle sets such as the high-resolution region of a zoom simulation.")
    end

    # integer cell coordinates
    icell = Matrix{Int32}(undef, 3, N)
    @inbounds for i = 1:N, k = 1:3
        icell[k, i] = floor(Int32, (pos[k, i] - lo[k]) / c)
    end
    n = vec(maximum(icell, dims=2)) .+ 1

    # sort particles by cell: columns (ix, iy) first, then iz within a column
    perm = collect(Int32, 1:N)
    for k = 3:-1:1
        perm = counting_sort(perm, view(icell, k, :), n[k])
    end
    x = pos[:, perm]

    # occupied cells, ordered by column and iz, and the first cell of every (ix, iy) column
    cell_first = Int32[]
    cell_iz = Int32[]
    col_first = zeros(Int32, n[1] * n[2] + 1)
    @inbounds for i = 1:N
        p, q = perm[i], perm[max(i - 1, 1)]
        if i == 1 || icell[1, p] != icell[1, q] || icell[2, p] != icell[2, q] || icell[3, p] != icell[3, q]
            push!(cell_first, i)
            push!(cell_iz, icell[3, p])
            col_first[icell[1, p]*n[2]+icell[2, p]+2] += 1
        end
    end
    push!(cell_first, N + 1)
    col_first[1] = 1
    for k = 2:length(col_first)
        col_first[k] += col_first[k-1]
    end

    Ncell = length(cell_iz)
    parent = collect(Int32, 1:Ncell)
    tree_size = ones(Int32, Ncell)

    # half of the 5×5 neighbour columns, the own column only links upwards in z
    col_offsets = ((0, 0), (0, 1), (1, -1), (1, 0), (1, 1),
                   (0, 2), (1, -2), (1, 2), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2))

    @inbounds for ix = 0:n[1]-1, iy = 0:n[2]-1

        col_a = ix * n[2] + iy + 1
        first_a, last_a = col_first[col_a], col_first[col_a+1] - 1
        first_a > last_a && continue

        for (dx, dy) ∈ col_offsets

            jx, jy = ix + dx, iy + dy
            (jx < n[1] && 0 <= jy < n[2]) || continue

            col_b = jx * n[2] + jy + 1
            first_b, last_b = col_first[col_b], col_first[col_b+1] - 1
            first_b > last_b && continue

            # walk both columns upwards in z
            p = first_b
            for a = first_a:last_a
                iz = cell_iz[a]
                while p <= last_b && cell_iz[p] < iz - 2
                    p += 1
                end
                b = p
                while b <= last_b && cell_iz[b] <= iz + 2
                    if dx != 0 || dy != 0 || cell_iz[b] > iz
                        link_cells!(parent, tree_size, a, b, x, cell_first, l2)
                    end
                    b += 1
                end
            end
        end
    end

    label = Vector{Int32}(undef, N)
    @inbounds for a = 1:Ncell
        r = find_root!(parent, a)
        for i = cell_first[a]:cell_first[a+1]-1
            label[perm[i]] = r
        end
    end

    return label
end

"""
    most_massive_fof_group(pos::AbstractMatrix{<:Real}, mass::AbstractVector{<:Real}, linking_length::Real)

Particle indices and total mass of the most massive friends-of-friends group of the
particles at `pos` with masses `mass`.
"""
function most_massive_fof_group(pos::AbstractMatrix{<:Real}, mass::AbstractVector{<:Real}, linking_length::Real)

    label = fof_labels(pos, linking_length)

    group_mass = zeros(length(label))
    @inbounds for i ∈ eachindex(label)
        group_mass[label[i]] += mass[i]
    end
    M, g = findmax(group_mass)

    return findall(==(g), label), M
end

"""
    shrinking_sphere_center(pos::AbstractMatrix{<:Real}, mass::AbstractVector{<:Real};
                            shrink::Real=0.975, n_min::Integer=1000)

Centre of a particle distribution after Power et al. (2003): the centre of mass of the
particles in a sphere that shrinks by the factor `shrink` per iteration around the current
centre. The iteration stops before the sphere holds fewer than 1% of the particles, but at
least 20 and at most `n_min` particles.
"""
function shrinking_sphere_center(pos::AbstractMatrix{<:Real}, mass::AbstractVector{<:Real};
                                 shrink::Real=0.975, n_min::Integer=1000)

    x = Float64.(pos)
    m = Float64.(mass)

    # the first n entries of idx are the particles in the sphere
    idx = collect(axes(x, 2))
    n = length(idx)
    iszero(n) && error("No particles to find the centre of.")
    n_stop = min(n_min, max(n ÷ 100, 20))

    center = shrinking_sphere_com(x, m, idx, n)
    r2 = 0.0
    @inbounds for i ∈ idx
        r2 = max(r2, (x[1, i] - center[1])^2 + (x[2, i] - center[2])^2 + (x[3, i] - center[3])^2)
    end

    # coincident particles end the iteration once the radius falls below the precision
    r2_end = eps()^2 * r2

    while r2 > r2_end
        r2 *= shrink^2

        # move the particles inside the smaller sphere to the front
        n_in = 0
        @inbounds for k = 1:n
            i = idx[k]
            if (x[1, i] - center[1])^2 + (x[2, i] - center[2])^2 + (x[3, i] - center[3])^2 < r2
                n_in += 1
                idx[k], idx[n_in] = idx[n_in], i
            end
        end
        n_in < n_stop && break

        n = n_in
        center = shrinking_sphere_com(x, m, idx, n)
    end

    return center
end

"""
    shrinking_sphere_com(x, m, idx, n)

Centre of mass of the particles `idx[1:n]`.
"""
function shrinking_sphere_com(x, m, idx, n)
    com = zeros(3)
    M = 0.0
    @inbounds for k = 1:n
        i = idx[k]
        com[1] += m[i] * x[1, i]
        com[2] += m[i] * x[2, i]
        com[3] += m[i] * x[3, i]
        M += m[i]
    end
    return com ./ M
end

"""
    fof_linking_length(snap_base::String, mass::AbstractVector{<:Real}; b::Real=0.2, omega_b::Real=0.04)

Linking length `b (m̄ / ρ̄_dm)^(1/3)` for particles with the masses `mass`, as in the
OpenGadget3 FoF. The mean matter density of the periodic box is the total particle mass over
the box volume, so no unit system is needed. Dark matter takes the share
`(Ω_0 - Ω_b) / Ω_0` of it, or all of it for a snapshot without gas.
"""
function fof_linking_length(snap_base::String, mass::AbstractVector{<:Real}; b::Real=0.2, omega_b::Real=0.04)

    h = read_header(snap_base)

    if iszero(h.boxsize)
        error("Snapshot has no periodic box to define the mean density, please provide `linking_length`.")
    end

    M = 0.0
    for parttype = 0:5
        iszero(get_total_particles(h, parttype)) && continue
        M += sum(Float64, read_block(snap_base, "MASS"; parttype))
    end

    Ω_dm = iszero(get_total_particles(h, 0)) ? h.omega_0 : h.omega_0 - omega_b
    ρ_dm = M / h.boxsize^3 * Ω_dm / h.omega_0

    return b * cbrt(sum(Float64, mass) / length(mass) / ρ_dm)
end

"""
    find_main_halo(snap_base::String; b::Real=0.2, omega_b::Real=0.04,
                   linking_length::Union{Real,Nothing}=nothing,
                   center::Symbol=:potential, verbose::Bool=true)

Centre of the most massive halo in the high-resolution region of a zoom simulation, in code
units. The halo is the most massive friends-of-friends group of the high-resolution dark
matter (particle type 1).

Only the high-resolution dark matter is linked. Gas particles that left the high-resolution
region, or particles next to massive low-resolution particles, can have a lower potential
than the halo centre, so the global potential minimum is not a reliable centre.

# Keyword Arguments
- `b`: linking length in units of the mean particle separation of type 1 at the mean dark
  matter density of the box (see [`GadgetIO.fof_linking_length`](@ref)).
- `omega_b`: baryon density parameter of the simulation, which is not stored in the header.
  Only used if the snapshot contains gas.
- `linking_length`: linking length in code units, replaces the one from `b`.
- `center`: definition of the centre
    - `:potential`: position of the group member with the lowest potential (needs the `POT` block)
    - `:shrinking_sphere`: shrinking-sphere centre of the group members (see [`GadgetIO.shrinking_sphere_center`](@ref))
- `verbose`: print the linking length and the properties of the group.

Positions are unwrapped around a high-resolution particle, so the centre is the periodic image
next to the high-resolution region, also if the halo lies across the border of the box.

# Examples
```julia
center = find_main_halo(snap_base)

# without a POT block in the snapshot
center = find_main_halo(snap_base, center=:shrinking_sphere)
```
"""
function find_main_halo(snap_base::String; b::Real=0.2, omega_b::Real=0.04,
                        linking_length::Union{Real,Nothing}=nothing,
                        center::Symbol=:potential, verbose::Bool=true)

    if center ∉ (:potential, :shrinking_sphere)
        error("center must be :potential or :shrinking_sphere, got :$center")
    end
    if center == :potential && !block_present(select_file(snap_base, 0), "POT")
        error("Block POT not present, use center=:shrinking_sphere for snapshots without potential.")
    end

    h = read_header(snap_base)

    pos = read_block(snap_base, "POS", parttype=1)
    mass = read_block(snap_base, "MASS", parttype=1)

    # unwrap the periodic box around the first particle
    if !iszero(h.boxsize)
        pos .= shift_across_box_border.(pos, pos[:, 1], h.boxsize, 1 // 2 * h.boxsize)
    end

    if isnothing(linking_length)
        linking_length = fof_linking_length(snap_base, mass; b, omega_b)
    end
    verbose && @info @sprintf("Linking length: %.4f", linking_length)

    members, M = most_massive_fof_group(pos, mass, linking_length)

    if center == :potential
        pot = read_block(snap_base, "POT", parttype=1)
        halo_center = Float64.(pos[:, members[argmin(pot[members])]])
    else
        halo_center = shrinking_sphere_center(pos[:, members], mass[members])
    end

    verbose && @info @sprintf("Main halo: %d particles, dark matter mass %.4e, centre (%.3f, %.3f, %.3f)",
                              length(members), M, halo_center...)

    return halo_center
end
