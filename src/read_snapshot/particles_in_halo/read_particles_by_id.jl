using Dates

"""
    filter_by_ids( selected_ids::Array{<:Integer}, ids::Array{<:Integer})

Checks for matching IDs.
"""
function filter_by_ids( selected_ids::Array{<:Integer}, ids::Array{<:Integer})

    # look for the array positions where the IDs match
    matched = sort( get_index_list( selected_ids, ids ) )

    return matched
end


"""
    read_filtered(snap_file, blockname, parttype, block_position, matched)

Helper function to return the correct filtered arrays.
"""
function read_filtered(snap_file, blockname, parttype, block_position, matched)
    data = read_block(snap_file, blockname, parttype=parttype, block_position=block_position)
    if size(data,2) > 1
        return data[:,matched]
    else
        return data[matched]
    end
end

"""
    read_particles_by_id_single_file(snap_file::String, halo_ids::Array{<:Integer}, 
                                          blocks::Array{String}, parttype::Integer)

Helper function to read the matching IDs of particles in one whole file.
"""
function read_particles_by_id_single_file(snap_file::String, halo_ids::Array{<:Integer}, 
                                          blocks::Array{String}, parttype::Integer; 
                                          verbose::Bool=true)

    # read ids in file
    ids = read_snap(snap_file, "ID", parttype) 

    if verbose
        @info "Matching IDs..."
        t1 = Dates.now()
    end

    # check which IDs match
    matched = filter_by_ids( halo_ids, ids)

    if verbose
        t2 = Dates.now()
        @info "Found $(size(matched,1)) matches. Took: $(t2 - t1)"
    end

    # store block positions for faster read-in
    block_positions = get_block_positions(snap_file)

    # read the data blocks whole, but only store relevant entries
    return Dict(blocks[i] => read_filtered(snap_file, blocks[i], parttype, block_positions[blocks[i]], matched) 
                for i = 1:size(blocks,1))

end

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
    blocks, no_mass_block = check_blocks(filename, blocks)

    # extend the list of blocks to read by ID block 
    blocks = [ blocks ; "ID" ]
    
    # remove duplicate blocks
    unique!(blocks)

    # for a given halo position and search radius we can use `read_particles_in_volume`
    if ( !isnothing(pos0) && r0 != 0.0)
        
        # check if .key files are present
        if isfile(filename * ".key")

            # read all particles in the defined volume
            data = read_particles_in_volume(snap_base, blocks, pos0, r0,
                                            parttype=parttype, verbose=verbose,
                                            use_keys=use_keys)

            if verbose
                println()
                @info "Matching IDs..."
                t1 = Dates.now()
            end

            # find matching entries
            matched = filter_by_ids( selected_ids, data["ID"])

            if verbose
                t2 = Dates.now()
                @info "Found $(size(matched,1)) matches. Took: $(t2 - t1)"
            end

            # construct a new dict only containing the matching particles
            return construct_matched_dict(data, blocks, matched)

        else # if there are no .key files we need to read the while snapshot
            return read_particles_by_id_single_file(filename, selected_ids, blocks, parttype, verbose=verbose)
        end

    else # otherwise we need to search whole files
       
        # single file read
        if isfile(snap_base)

            return read_particles_by_id_single_file(snap_base, selected_ids, blocks, parttype, verbose=verbose)
            
        else # multi-file brute-force read -> slow!

            # total number of particles to read
            N_to_read = size(selected_ids,1)

            # number of particles read so far 
            N_read = 0

            # read the info block
            snap_info = read_info(filename)

            # prepare data dict for storage
            data = allocate_data_dict(blocks, N_to_read, snap_info, no_mass_block)

            filename = snap_base * ".0"
            h = read_header(filename)

            @warn "Brute-force reading $(h.num_files) files! This may take a while!"

            # loop over all the files in the snap directory
            for i = 0:(h.num_files-1)

                # get current file name
                filename = select_file(snap_base, i)

                # read data from file
                data_file = read_particles_by_id_single_file(filename, selected_ids, blocks, parttype, verbose=verbose)

                N_this_file = size(data_file["ID"],1)

                # write into master dict
                for block ∈ blocks
                    dim   = snap_info[getfield.(snap_info, :block_name) .== block][1].n_dim
                    if dim == 1
                        data[block][N_read+1:N_read+N_this_file]    = data_file[block]
                    else
                        data[block][:, N_read+1:N_read+N_this_file] = data_file[block]
                    end
                end # blocks

                # update number of read particles
                N_read += N_this_file

            end

            # reduce array size
            for block ∈ blocks
                # eltype needed for type stable reshape function call
                dim   = eltype(N_read)(snap_info[getfield.(snap_info, :block_name) .== block][1].n_dim)
                if dim == 1
                    resize!(data[block], N_read)
                else
                    data[block] = reshape(resize!(vec(data[block]), dim*N_read), dim, N_read)
                end
            end # blocks

            return data

        end
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
