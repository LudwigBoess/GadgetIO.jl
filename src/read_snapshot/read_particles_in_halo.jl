using Dates

"""
    select_subfile(sub_base::String, filenum::Integer)

Checks if the requested files exists and returns the path.
"""
function select_subfile(sub_base::String, filenum::Integer)

    if !isfile(sub_base)
        return sub_base * ".$filenum"
    else
        return sub_base
    end
end

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
    read_halo_pids(sub_base::String, offset::Integer, N_ids::Integer)

Reads the IDs of the particles that are contained in the halo.
"""
function read_halo_pids(sub_base::String, offset::Integer, N_ids::Integer)

    sub_file = select_subfile(sub_base, 0)

    sub_header = read_subfind_header(sub_file)

    # if there are fewer IDs in the file than in the halo they are distributed over multiple files
    if (sub_header.nfof - offset) < N_ids
        
        # allocate array for IDs
        ids = Array{<:Integer}(undef, N_ids)

        ids_read = 0
        files_read = 0

        # the remainig offset
        offset_remaining = offset

        # loop over sub files untiul all IDs are read
        while ids_read < N_ids

            sub_file = select_subfile(sub_base, files_read)

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
    read_particles_in_halo(snap_base::String, blocks::Array{String}, 
                           sub_base::String, halo::SubfindID, halo_type::Integer)

Reads all particles of type `parttype` that are contained in a halo defined by its `SubfindID`.
Returns a `Dict` with each of the `blocks` as entries.
"""
function read_particles_in_halo(snap_base::String, blocks::Array{String},
                                sub_base::String, halo::SubfindID, halo_type::Integer;
                                parttype::Integer=0, verbose::Bool=true)

    # select subfind file to read
    sub_file = select_subfile(sub_base, halo.file)

    # select block with number of particles in fof for
    if halo_type == 1
        # halos
        len_block = "GLEN"
        off_block = "GOFF"
        pos_block = "GPOS"
        rad_block = "R200"
        rad_scale = 1.0
    elseif halo_type == 2
        # subhalos
        len_block = "SLEN"
        off_block = "SOFF"
        pos_block = "SPOS"
        rad_block = "RHMS"
        rad_scale = 5.0
    end

    if verbose
        @info "Reading number of particles in halo..."
        t1 = Dates.now()
    end
    
    # read number of IDs in the halo
    N_ids = read_subfind(sub_file, len_block)[halo.id]

    if verbose
        t2 = Dates.now()
        @info "N_particles to read: $N_ids. Took: $(t2 - t1)"
    end
    # read offset in PID array
    offset = read_subfind(sub_file, off_block)[halo.id]

    if verbose
        @info "Reading IDs in halo..."
        t1 = Dates.now()
    end
    # read all IDs of the particles contained in a halo
    halo_ids = read_halo_pids(sub_base, offset, N_ids)

    if verbose
        t2 = Dates.now()
        @info "IDs read. Took: $(t2 - t1)"

        @info "Reading Data..."
        t1 = Dates.now()
    end

    if !isfile(snap_base)
        # position of halo
        halo_pos = read_subfind(sub_file, pos_block)[halo.id,:]

        # initial search radius for read-in
        initial_radius = rad_scale * read_subfind(sub_file, rad_block)[halo.id]

        blocks = [ blocks ; "ID" ]

        data = read_particles_in_volume(snap_base, blocks, halo_pos, initial_radius,
                                        parttype=parttype, verbose=verbose)

        if verbose
            t2 = Dates.now()
            println()
            @info "Data read. Took: $(t2 - t1)"

            @info "Matching IDs..."
            t1 = Dates.now()
        end

        matched = filter_by_ids( halo_ids, data["ID"])

        if verbose
            t2 = Dates.now()
            @info "IDs found. Took: $(t2 - t1)"
            println()
            @info "Done!"
        end


        return Dict(blocks[i] => data[blocks[i]][matched,:] for i = 1:length(blocks))

    else # only read one file

        # read all blocks into a dictionary
        data = Dict(blocks[i] => read_snap(snap_base, blocks[i], parttype) for i = 1:length(blocks))

        if verbose
            t2 = Dates.now()
            println()
            @info "Data read. Took: $(t2 - t1)"

            @info "Matching IDs..."
            t1 = Dates.now()
        end

        # check which IDs match
        matched = filter_by_ids( halo_ids, data["ID"])

        if verbose
            t2 = Dates.now()
            @info "IDs found. Took: $(t2 - t1)"
            println()
            @info "Done!"
        end

        return Dict(blocks[i] => data[blocks[i]][matched,:] for i = 1:length(blocks))
    end

end