"""
    read_particles_in_geometry( filename::String, blocks::Vector{String},
                                geometry::AbstractGadgetGeometry;
                                parttype::Integer=0, verbose::Bool=true,
                                use_keys::Bool=true)

Reads all particles within a space defined by an `AbstractGeometry` struct for a given particle type. 
Returns a dictionary with all requested blocks.
"""
function read_particles_in_geometry(filename::String, blocks::Vector{String},
                                    geometry::AbstractGadgetGeometry;
                                    parttype::Integer=0, verbose::Bool=true,
                                    use_keys::Bool=true)

    # make sure position is read - needed for determining if in geometry or not
    if !("POS" ∈ blocks)
        @info "POS not found in blocks, adding POS..."
        append!(blocks, "POS")
    end


    # get the corners of a box which contains the geometry
    corner_lowerleft, corner_upperright = get_geometry_box_corners(geometry)

    # read the box
    d = read_particles_in_box(filename, blocks, corner_lowerleft, corner_upperright; parttype, verbose, use_keys)

    # determine particles within geometry
    mask = get_geometry_mask(geometry, d["POS"])

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