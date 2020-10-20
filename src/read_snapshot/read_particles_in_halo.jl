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
        @info "Found $(length(matched)) matches. Took: $(t2 - t1)"
    end

    # read the data blocks whole, but only store relevant entries
    return Dict(blocks[i] => read_snap(snap_file, blocks[i], parttype)[matched,:] for i = 1:length(blocks))

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
                @info "Found $(length(matched)) matches. Took: $(t2 - t1)"
            end
            
            return Dict(blocks[i] => data[blocks[i]][matched,:] for i = 1:length(blocks))

        else # if there are no .key files we need to read the while snapshot
            return read_particles_by_id_single_file(snap_file, selected_ids, blocks, parttype, verbose=verbose)
        end

    else # otherwise we need to search whole files
       
        # single file read
        if isfile(snap_base)

            return read_particles_by_id_single_file(snap_base, selected_ids, blocks, parttype, verbose=verbose)
            
        else # multi-file brute-force read -> slow!

            # total number of particles to read
            N_to_read = length(ids)

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

                N_this_file = length(data_file["ID"][:,1])

                # write into master dict
                for block in blocks
                    data[block][N_read+1:N_read+N_this_file,:] = data_file[block]
                end # blocks

                # update number of read particles
                N_read += N_this_file

            end

            # reduce array size
            for block in blocks
                if length(data[block][1,:]) == 1
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


"""
    read_pids(sub_base::String, offset::Integer, N_ids::Integer)

Reads the `PID` block in the subfind output.
"""
function read_pids(sub_base::String, N_ids::Integer, offset::Integer)

    # choose correct file
    sub_file = select_file(sub_base, 0)

    # read the header
    sub_header = read_subfind_header(sub_file)

    # read info for PID datatype
    info = read_info(sub_file)
    pid_info = info[getfield.(info, :block_name) .== "PID"][1]

    # if there are fewer IDs in the file than in the halo they are distributed over multiple files
    if (sub_header.nfof - offset) < N_ids
        
        # allocate array for IDs
        ids = Array{pid_info.data_type}(undef, N_ids)

        ids_read = 0
        files_read = 0

        # the remainig offset
        offset_remaining = offset

        # loop over sub files untiul all IDs are read
        while ids_read < N_ids

            sub_file = select_file(sub_base, files_read)

            if !isfile(sub_file)
                error("$sub_file does not exist!")
            end

            # read the header
            sub_header = read_subfind_header(sub_file)

            # if the offset is larger than the number of IDs in this file
            if (sub_header.nfof - offset_remaining) < 0.0
                # read next file
                files_read += 1
                # subtract IDs in this file from offset
                offset_remaining -= sub_header.nfof
                continue
            else
                # find starting position for read-in
                start_pos = offset_remaining + 1
                
                # number of particles remaining to be read
                if ( sub_header.nfof - offset_remaining < N_ids - ids_read )
                    n_to_read = sub_header.nfof - offset_remaining
                else
                    n_to_read = N_ids - ids_read
                end

                ids[ids_read+1:ids_read+n_to_read] = read_subfind(sub_file, "PID")[start_pos:start_pos+n_to_read-1]

                # update number of IDs alrady read
                ids_read += n_to_read

                # increase for next file read
                files_read += 1
            end
            
        end # while

        return ids
    else # they can be read from one file
        # position of first ID in array
        start_pos = offset + 1
        return read_subfind(sub_file, "PID")[start_pos:start_pos+N_ids-1]
    end
end


"""
    read_ids_in_halo( sub_base::String, halo::HaloID; 
                      halo_type::Integer=1, verbose::Bool=true)

Reads the IDs of all particles contained in a `halo`.
"""
function read_ids_in_halo( sub_base::String, halo::HaloID; 
                           halo_type::Integer=1, verbose::Bool=true)


    # select subfind file to read
    sub_file = select_file(sub_base, halo.file)

    # select block with number of particles in fof for
    if halo_type == 1
        # halos
        len_block = "GLEN"
        off_block = "GOFF"

    elseif halo_type == 2
        # subhalos
        len_block = "SLEN"
        off_block = "SOFF"
    end

    if verbose
        @info "Reading number of particles in halo..."
        t1 = Dates.now()
    end
    
    # read number of IDs in the halo
    N_ids = Int64(read_subfind(sub_file, len_block)[halo.id])

    if verbose
        t2 = Dates.now()
        @info "N_particles to read: $N_ids. Took: $(t2 - t1)"
    end
    # read offset in PID array
    offset = Int64(read_subfind(sub_file, off_block)[halo.id])

    if verbose
        @info "Reading IDs in halo..."
        t1 = Dates.now()
    end
    # read all IDs of the particles contained in a halo
    halo_ids = read_pids(sub_file, N_ids, offset)

    if verbose
        t2 = Dates.now()
        @info "IDs read. Took: $(t2 - t1)"
    end

    return sort(halo_ids)
end


"""
    read_particles_in_halo(snap_base::String, blocks::Array{String},
                                sub_base::String, halo::HaloID; 
                                rad_scale::Real=1.0, halo_type::Integer=1,
                                parttype::Integer=0, verbose::Bool=true)

Reads all particles of type `parttype` that are contained in a halo defined by its `HaloID`.
Returns a `Dict` with each of the `blocks` as entries.
"""
function read_particles_in_halo(snap_base::String, blocks::Array{String},
                                sub_base::String, halo::HaloID; 
                                rad_scale::Real=1.0, halo_type::Integer=1,
                                parttype::Integer=0, verbose::Bool=true)

    # select subfind file to read
    sub_file = select_file(sub_base, halo.file)

    # select block with number of particles in fof for
    if halo_type == 1
        # halos
        pos_block = "GPOS"
        rad_block = "R200"
    elseif halo_type == 2
        # subhalos
        pos_block = "SPOS"
        rad_block = "RHMS"
    end

    if verbose
        @info "Reading IDs in halo..."
        t1 = Dates.now()
    end
    # read all IDs of the particles contained in a halo
    halo_ids = read_ids_in_halo(sub_base, halo, 
                                halo_type=halo_type, verbose=verbose)

    if verbose
        t2 = Dates.now()
        @info "IDs read. Took: $(t2 - t1)"

        @info "Reading Data..."
        t1 = Dates.now()
    end

    # position of halo
    halo_pos = read_subfind(sub_file, pos_block)[halo.id,:]

    # initial search radius for read-in
    initial_radius = rad_scale * read_subfind(sub_file, rad_block)[halo.id]

    # read ids in the halo 
    data = read_particles_by_id(snap_base, halo_ids, blocks,
                              parttype=parttype, verbose=verbose,
                              pos0=halo_pos,
                              r0=initial_radius)

    if verbose
        t2 = Dates.now()
        @info "Data read. Took: $(t2 - t1)"
        println()
        @info "Done!"
    end

    return data

end


"""
    read_particles_in_halo( snap_base::String, block::String,
                            sub_base::String, halo::HaloID; 
                            rad_scale::Real=1.0, halo_type::Integer=1,
                            parttype::Integer=0, verbose::Bool=true )

Reads all particles of type `parttype` that are contained in a halo defined by its `HaloID`.
Returns an `Array` with the requested `block`.
"""
function read_particles_in_halo(snap_base::String, block::String,
                                sub_base::String, halo::HaloID; 
                                rad_scale::Real=1.0, halo_type::Integer=1,
                                parttype::Integer=0, verbose::Bool=true)

    data = read_particles_in_halo(snap_base, [block],
                                  sub_base, halo, 
                                  rad_scale=rad_scale, halo_type=halo_type,
                                  parttype=parttype, verbose=verbose)

    return data[block]
end
