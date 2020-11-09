using Dates

"""
    filter_by_ids( selected_ids::Array{<:Integer}, ids::Array{<:Integer})

Checks for matching IDs.
"""
function filter_by_ids( selected_ids::Array{<:Integer}, ids::Array{<:Integer})

    # look for the array positions where the IDs match
    matched = sort( get_index_list_dict( selected_ids, ids ) )

    return matched
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
        @info "Found $(size(matched)[1]) matches. Took: $(t2 - t1)"
    end

    # read the data blocks whole, but only store relevant entries
    return Dict(blocks[i] => read_snap(snap_file, blocks[i], parttype)[matched,:] for i = 1:size(blocks)[1])

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
                              pos0::Array{<:Real}=[-1.234, -1.234, -1.234],
                              r0::Real=0.0,
                              use_keys::Bool=true)

    # sort the IDs if they are not already sorted
    if !issorted(selected_ids)
        sort!(selected_ids)
    end
    # try reading the first of the distributed snapshots
    filename = select_file(snap_base, 0)

    # check if blocks are present 
    blocks, no_mass_block = check_blocks(filename, blocks)

    # extend the list of blocks to read by ID block 
    blocks = [ blocks ; "ID" ]
    unique!(blocks)

    # for a given halo position and search radius we can use `read_particles_in_volume`
    if (pos0 != [-1.234, -1.234, -1.234] && r0 != 0.0)
        
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
                @info "Found $(size(matched)[1]) matches. Took: $(t2 - t1)"
            end
            
            return Dict(blocks[i] => data[blocks[i]][matched,:] for i = 1:size(blocks)[1])

        else # if there are no .key files we need to read the while snapshot
            return read_particles_by_id_single_file(filename, selected_ids, blocks, parttype, verbose=verbose)
        end

    else # otherwise we need to search whole files
       
        # single file read
        if isfile(snap_base)

            return read_particles_by_id_single_file(snap_base, selected_ids, blocks, parttype, verbose=verbose)
            
        else # multi-file brute-force read -> slow!

            # total number of particles to read
            N_to_read = size(ids)[1]

            # number of particles read so far 
            N_read = 0

            # read the info block
            snap_info = read_info(filename)

            # prepare data dict for storage
            data = allocate_data_dict(blocks, N_to_read, snap_info, no_mass_block)

            filename = snap_base * ".0"
            h = head_to_obj(filename)
            nfiles = h.num_files

            @warn "Brute-force reading $numfiles files! This may take a while!"

            # loop over all the files in the snap directory
            for i = 0:nfiles

                # get current file name
                filename = select_file(snap_base, i)

                # read data from file
                data_file = read_particles_by_id_single_file(filename, selected_ids, blocks, parttype, verbose=verbose)

                N_this_file = size(data_file["ID"])[1]

                # write into master dict
                for block in blocks
                    data[block][N_read+1:N_read+N_this_file,:] = data_file[block]
                end # blocks

                # update number of read particles
                N_read += N_this_file

            end

            # reduce array size
            for block in blocks
                if size(data[block])[2] == 1
                    resize!(data[block], N_read)
                else
                    data[block] = reshape(resize!(vec(data[block]),3*N_read),N_read,3)
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
                              pos0::Array{<:Real}=[-1.234, -1.234, -1.234],
                              r0::Real=0.0)

   data = read_particles_by_id(snap_base, selected_ids, 
                              [block],
                              parttype=parttype, verbose=verbose,
                              pos0=pos0, r0=r0)

    return data[block]
end