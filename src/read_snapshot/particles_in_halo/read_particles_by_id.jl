using Dates


"""
    construct_matched_dict(data_in::Dict{String, Union{Vector{T}, Array{T,2}}}, 
                                blocks::Array{String,1}, matched::Array{<:Integer,1}) where T

Write all matching particles to a new dictionary.
"""
function construct_matched_dict(data_in::Dict{String, VecOrMat{T} where T}, 
                                blocks::Array{String,1}, matched::Array{<:Integer,1})

    dict_out = Dict{String, VecOrMat{T} where T}()

    for block ∈ blocks
        dim   = size(data_in[block], 2)
        if dim == 1
            dict_out[block] = data_in[block][matched]
        else
            dict_out[block] = data_in[block][:, matched]
        end
    end
    return dict_out
end


"""
    read_particles_by_id(snap_base::String, ids::Array{<:Integer}, 
                         blocks::Array{String}; 
                         parttype::Integer=0, verbose::Bool=true,
                         pos0::Array{<:Real}=[-1.234, -1.234, -1.234],
                         r0::Real=0.0)

Reads particles filtered by the provided IDs. Returns all requested blocks as entries in a `Dict`.
"""
function read_particles_by_id(snap_base::String, selected_ids::Array{<:Integer}, 
                              blocks::Array{String}; 
                              parttype::Integer=0, verbose::Bool=true,
                              pos0::Union{Array{<:Real},Nothing}=nothing,
                              r0::Real=0.0,
                              use_keys::Bool=true)

    # try reading the first of the distributed snapshots
    filename = select_file(snap_base, 0)

    # check if blocks are present 
    blocks, no_mass_block = check_blocks(filename, blocks, parttype)

    # extend the list of blocks to read by ID block 
    blocks = [ blocks ; "ID" ]
    
    # remove duplicate blocks
    unique!(blocks)

    # for a given halo position and search radius we can use `read_particles_in_volume`
    if ( !isnothing(pos0) && r0 != 0.0)
        
        # read all particles in the defined volume
        data = read_particles_in_volume(snap_base, blocks, pos0, r0;
                                        parttype, verbose, use_keys)

        if verbose
            println()
            @info "Matching IDs..."
            t1 = Dates.now()
        end

        # find matching entries
        matched = get_index_list( selected_ids, data["ID"] )

        if verbose
            t2 = Dates.now()
            @info "Found $(size(matched,1)) matches. Took: $(t2 - t1)"
        end

        # construct a new dict only containing the matching particles
        return construct_matched_dict(data, blocks, matched)

    else # otherwise we need to search whole files
        
        filter_function(filename) = filter_ids(filename, selected_ids, parttype)
        return read_blocks_filtered(snap_base, blocks; 
                                          filter_function,
                                          parttype, verbose)
        
    end
end

"""
    read_particles_by_id(snap_base::String, ids::Array{<:Integer}, 
                         block::String; 
                         parttype::Integer=0, verbose::Bool=true,
                         pos0::Array{<:Real}=[-1.234, -1.234, -1.234],
                         r0::Real=0.0)

Reads particles filtered by the provided IDs. Returns the requested block as an `Array`.
"""
function read_particles_by_id(snap_base::String, selected_ids::Array{<:Integer}, 
                              block::String; 
                              parttype::Integer=0, verbose::Bool=true,
                              pos0::Union{Array{<:Real},Nothing}=nothing,
                              r0::Real=0.0,
                              use_keys::Bool=true)

   data = read_particles_by_id(snap_base, selected_ids, 
                              [block],
                              parttype=parttype, verbose=verbose,
                              pos0=pos0, r0=r0, use_keys=use_keys)

    return data[block]
end
