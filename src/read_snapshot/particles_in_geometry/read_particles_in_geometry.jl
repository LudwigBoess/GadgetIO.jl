"""
    read_particles_in_geometry( filename::String, blocks::Vector{String},
                                geometry::AbstractGadgetGeometry;
                                parttype::Integer=0, verbose::Bool=true,
                                use_keys::Bool=true, shift_across_box_border::Bool=true)

Reads all particles within a space defined by an `AbstractGeometry` struct for a given particle type. 
Returns a dictionary with all requested blocks.
If `origin_to_geometry` is `true`, the origin is moved to the geometry's defined center, respecting any
periodicity of the box.
"""
function read_particles_in_geometry(filename::String, blocks::Vector{String},
                                    geometry::AbstractGadgetGeometry;
                                    parttype::Integer=0, verbose::Bool=true,
                                    use_keys::Bool=true, shift_across_box_border::Bool=true)

    h = read_header(filename * ".0")

    # make sure position is read - needed for determining if in geometry or not
    if !("POS" ∈ blocks)
        @info "POS not found in blocks, adding POS..."
        append!(blocks, "POS")
    end


    # get the corners of a box which contains the geometry
    corner_lowerleft, corner_upperright = get_geometry_box_corners(geometry)

    # read the box
    d = read_particles_in_box(filename, blocks, corner_lowerleft, corner_upperright; parttype, verbose, use_keys)

    # determine particles within geometry, shifting the particles across periodic boundaries
    r₀ = get_geometry_center(geometry)
    if shift_across_box_border
        d["POS"] .= GadgetIO.shift_across_box_border.(d["POS"], r₀, h.boxsize, 1 // 2 * h.boxsize)
        mask = get_geometry_mask(geometry, d["POS"])
    else
        pos = GadgetIO.shift_across_box_border.(d["POS"], r₀, h.boxsize, 1 // 2 * h.boxsize)
        mask = get_geometry_mask(geometry, pos)
    end

    # iterate blocks of dict to remove filtered out particles
    for block ∈ keys(d)
        # get all dimensions but the last (the last is masked)
        colons = repeat([:], ndims(d[block]) - 1)

        d[block] = d[block][colons..., mask]        
    end

    return d
end



"""
    read_particles_in_geometry( filename::String, block::String,
                                geometry::AbstractGadgetGeometry;
                                parttype::Integer=0, verbose::Bool=true,
                                use_keys::Bool=true)

Reads all particles within a space defined by an `AbstractGeometry` struct for a given particle type. 
Returns a dictionary with the requested block.
"""
function read_particles_in_geometry(filename::String, block::String,
                                    geometry::AbstractGadgetGeometry;
                                    parttype::Integer=0, verbose::Bool=true,
                                    use_keys::Bool=true)

    return read_particles_in_geometry(filename, [block], geometry; parttype, verbose, use_keys)
end


""" 
    shift_across_box_border(x::Real, x_halo::Real, boxsize::Real, boxsize_half::Real)

Shift coordinate `x` across the box border if the zero coordinate `x₀` is on the other side.
"""
function shift_across_box_border(x::Real, x₀::Real, boxsize::Real, boxsize_half::Real)
    if x - x₀ > boxsize_half
        return x - boxsize
    elseif x₀ - x > boxsize_half
        return x + boxsize
    end
    return x
end
