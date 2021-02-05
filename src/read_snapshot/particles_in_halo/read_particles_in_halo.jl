using Dates

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

        # loop over sub files until all IDs are read
        while ids_read < N_ids

            sub_file = select_file(sub_base, files_read)

            if !isfile(sub_file)
                error("$sub_file does not exist!")
            end

            # read the header
            sub_header = read_subfind_header(sub_file)

            # if the offset is larger than the number of IDs in this file
            if offset_remaining >= sub_header.nfof
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
    halo_ids = read_pids(sub_base, N_ids, offset)

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
                                parttype::Integer=0, verbose::Bool=true,
                                use_keys::Bool=true)

    # select subfind file to read
    sub_file = select_file(sub_base, halo.file)

    # select block with number of particles in fof for
    if halo_type == 1
        # halos
        pos_block = "GPOS"
        rad_block = "R200"
        if !block_present(sub_file, rad_block)
            rad_block = "RMEA"

            if !block_present(sub_file, rad_block)
                error("Neither R200 nor RMEA present!")
            end
        end

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
    halo_pos = read_subfind(sub_file, pos_block)[:,halo.id]

    # initial search radius for read-in
    initial_radius = rad_scale * read_subfind(sub_file, rad_block)[halo.id]

    # read ids in the halo 
    data = read_particles_by_id(snap_base, halo_ids, blocks,
                              parttype=parttype, verbose=verbose,
                              pos0=halo_pos,
                              r0=initial_radius,
                              use_keys=use_keys)

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
                                parttype::Integer=0, verbose::Bool=true,
                                use_keys::Bool=true)

    data = read_particles_in_halo(snap_base, [block],
                                  sub_base, halo, 
                                  rad_scale=rad_scale, halo_type=halo_type,
                                  parttype=parttype, verbose=verbose,
                                  use_keys=use_keys)

    return data[block]
end
