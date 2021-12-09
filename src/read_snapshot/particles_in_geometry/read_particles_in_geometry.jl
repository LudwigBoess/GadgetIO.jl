"""
    read_particles_in_geometry( filename::String, blocks::Vector{String},
                                geometry::AbstractGadgetGeometry;
                                parttype::Integer=0, verbose::Bool=true,
                                use_keys::Bool=true, do_shift_across_box_border::Bool=true)

Reads all particles within a space defined by an `AbstractGeometry` struct for a given particle type. 
Returns a dictionary with all requested blocks.
If `shift_across_box_border` is `true`, the particles are moved beyond the borders of a periodic box,
if `false`, the periodicity of the box is still accounted for, but the particles' positions are not
shifted.
"""
function read_particles_in_geometry(filename::String, blocks::Vector{String},
                                    geometry::AbstractGadgetGeometry;
                                    parttype::Integer=0, verbose::Bool=true,
                                    use_keys::Bool=true, do_shift_across_box_border::Bool=true)

    
    h = read_header( select_file(filename, 0) )

    # make sure position is read - needed for determining if in geometry or not
    if "POS" ∉ blocks
        @info "POS not found in blocks, adding POS..."
        push!(blocks, "POS")
    end


    # get the corners of a box which contains the geometry
    corner_lowerleft, corner_upperright = get_geometry_box_corners(geometry)

    # read the box
    d = read_particles_in_box(filename, blocks, corner_lowerleft, corner_upperright; parttype, verbose, use_keys)

    # determine particles within geometry, shifting the particles across periodic boundaries
    r₀ = get_geometry_center(geometry)
    if do_shift_across_box_border
        d["POS"] .= shift_across_box_border.(d["POS"], r₀, h.boxsize, 1 // 2 * h.boxsize)
        mask = get_geometry_mask(geometry, d["POS"])
    else
        pos = shift_across_box_border.(d["POS"], r₀, h.boxsize, 1 // 2 * h.boxsize)
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
                                use_keys::Bool=true, do_shift_across_box_border::Bool=true)

Reads all particles within a space defined by an `AbstractGeometry` struct for a given particle type. 
Returns a dictionary with the requested block.
"""
function read_particles_in_geometry(filename::String, block::String,
                                    geometry::AbstractGadgetGeometry;
                                    parttype::Integer=0, verbose::Bool=true,
                                    use_keys::Bool=true, do_shift_across_box_border::Bool=true)

    return read_particles_in_geometry(filename, [block], geometry; parttype, verbose, use_keys)
end
