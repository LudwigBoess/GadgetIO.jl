"""
    These functions read particles in a given range over multiple files, based
    on the relevant peano-hilbert keys.
    This code is based on read_particles_in_box.pro by Dr. habil. Klaus Dolag.
"""

using Dates
using Base.Threads



"""
    read_particles_in_box_peano(filename::String, blocks::Vector{String},
                                corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                                parttype::Integer=0, verbose::Bool=true)

Reads all particles within a box defined by a lower left and upper right corner
for a given particle type based on peano hilbert key reading. Returns a dictionary with all requested blocks.
"""
function read_particles_in_box_peano(snap_base::String, blocks::Vector{String},
                                    corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                                    parttype::Integer=0, verbose::Bool=true)

    # get read positions from PH keys
    read_positions = read_positions_from_PH_keys(snap_base, corner_lowerleft, corner_upperright;
                                                 parttype, verbose)

    # call read function
    return read_blocks_filtered(snap_base, blocks; read_positions, parttype, verbose)
end


"""
    read_particles_in_box(filename::String, blocks::Vector{String},
                        corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                        parttype::Integer=0, verbose::Bool=true,
                        use_keys::Bool=true)

Reads all particles within a box defined by a lower left and upper right corner
for a given particle type. Returns a dictionary with all requested blocks.
If `use_keys=true` it uses Peano-Hilbert keys to constrain the read-in, otherwise it uses a brute-force read-in with a filter function.
Peano-Hilbert key based read-in is significantly faster.
"""
function read_particles_in_box(filename::String, blocks::Vector{String},
                               corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                               parttype::Integer=0, verbose::Bool=true,
                               use_keys::Bool=true)

    if use_keys
        d = read_particles_in_box_peano(filename, blocks, corner_lowerleft, corner_upperright, 
                                        parttype=parttype, verbose=verbose)
    else
        if verbose
            @info "Brute-force read-in."
        end
        filter_function(snap_file) = filter_cube(snap_file, corner_lowerleft, corner_upperright, parttype=parttype)
        d = read_blocks_filtered(filename, blocks, filter_function = filter_function, parttype = parttype, verbose = verbose )
    end

    return d
end

"""
    read_particles_in_box(filename::String, block::String,
                          corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                          parttype::Integer=0, verbose::Bool=true,
                          use_keys::Bool=true)

Like `read_particles_in_box` but for a single block. Returns the block as an array.

See also: [`read_particles_in_box`](@ref)
"""
function read_particles_in_box(filename::String, block::String,
                               corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                               parttype::Integer=0, verbose::Bool=true,
                               use_keys::Bool=true)

    d = read_particles_in_box(filename, [block], corner_lowerleft, corner_upperright, 
                              parttype=parttype, verbose=verbose, use_keys=use_keys)

    return d[block]
end


"""
    read_particles_in_box(filename::String, blocks::Vector{String},
                          center_pos, radius;
                          parttype::Integer=0, verbose::Bool=true,
                          use_keys::Bool=true)

Reads all particles within a box encapsulating a volume defined by center position
and radius for a given particle type. Returns a dictionary with all requested blocks.

See also: [`read_particles_in_box`](@ref)
"""
function read_particles_in_volume(filename::String, blocks::Vector{String},
                                  center_pos::Array{<:Real}, radius::Real;
                                  parttype::Integer=0, verbose::Bool=true,
                                  use_keys::Bool=true)

    # calculate lower left and upper right corner
    x0 = center_pos .- radius
    x1 = center_pos .+ radius

    return read_particles_in_box(filename, blocks, x0, x1, parttype=parttype, verbose=verbose, use_keys=use_keys)
end

"""
    read_particles_in_box(filename::String, block::String,
                          center_pos, radius;
                          parttype::Integer=0, verbose::Bool=true)

Like `read_particles_in_volume` but for a single block. Returns the block as an array.

See also: [`read_particles_in_volume`](@ref)
"""
function read_particles_in_volume(filename::String, block::String,
                                  center_pos, radius;
                                  parttype::Integer=0, verbose::Bool=true,
                                  use_keys::Bool=true)

    d = read_particles_in_volume(filename, [block], center_pos, radius,
                                 parttype=parttype, verbose=verbose,
                                 use_keys=use_keys)

    return d[block]
end

