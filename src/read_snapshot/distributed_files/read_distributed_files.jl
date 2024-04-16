using ProgressMeter
using Base.Threads

"""
    filter_positions(snap_file::String, corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real}, 
                          parttype::Integer)

Reads positions from `snap_file` and returns the indices of particles contained in a box defined by `corner_lowerleft` and `corner_upperright`.
"""
function filter_positions(snap_file::String, corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real}, 
                          parttype::Integer)

    if HDF5.ishdf5(snap_file)
        pos = h5read(snap_file, "/PartType$parttype")["Coordinates"]
    else
        pos = read_block(snap_file, "POS"; parttype)
    end

    sel = findall( @views ( ( corner_lowerleft[1] .<= pos[1,:] .< corner_upperright[1] ) .&
                            ( corner_lowerleft[2] .<= pos[2,:] .< corner_upperright[2] ) .&
                            ( corner_lowerleft[3] .<= pos[3,:] .< corner_upperright[3] )) )
    return sel
end

"""
    filter_ids(snap_file::String, selected_ids::Vector{<:Integer})

Reads IDs from `snap_file` and returns the indices of particles matching the `selected_ids`
"""
function filter_ids(snap_file::String, selected_ids::Vector{<:Integer}, parttype::Integer)

    if HDF5.ishdf5(snap_file)
        IDs = h5read(snap_file, "/PartType$parttype")["ParticleIDs"]
    else
        IDs = read_block(snap_file, "ID"; parttype)
    end

    return get_index_list(selected_ids, IDs)
end


"""
    read_blocks_filtered( snap_base::String, blocks::Array{String};
                                filter_function::Union{Function, Nothing}=nothing, 
                                read_positions::Union{Dict, Nothing}=nothing, 
                                parttype::Integer=0, verbose::Bool=true )

Reads the specified blocks from all distributed files where particles pass the `filter_function`, or are given by a `Dict` of `read_positions`.
For `read_positions` please see [`find_read_positions`](@ref).
"""
function read_blocks_filtered(snap_base::String, blocks::Array{String};
                              filter_function::Union{Function, Nothing}=nothing, 
                              read_positions::Union{Dict, Nothing}=nothing, 
                              parttype::Integer=0, verbose::Bool=true )

    # select a dummy file
    snap_test = select_file(snap_base, 0)   

    # check if the file is hdf5
    if HDF5.ishdf5(snap_test)
        return read_blocks_filtered_hdf5(snap_base, blocks; filter_function, read_positions, parttype, verbose)
    else
        return read_blocks_filtered_format2(snap_base, blocks; filter_function, read_positions, parttype, verbose)
    end
end

"""
    read_blocks_filtered_format2( snap_base::String, blocks::Array{String};
                                filter_function::Union{Function, Nothing}=nothing, 
                                read_positions::Union{Dict, Nothing}=nothing, 
                                parttype::Integer=0, verbose::Bool=true )

Reads the specified blocks from all distributed files where particles pass the `filter_function`, or are given by a `Dict` of `read_positions`.
For `read_positions` please see [`find_read_positions`](@ref).
"""
function read_blocks_filtered_format2(snap_base::String, blocks::Array{String};
                                      filter_function::Union{Function, Nothing}=nothing, 
                                      read_positions::Union{Dict, Nothing}=nothing, 
                                      parttype::Integer=0, verbose::Bool=true )

    # default behaviour is to find the particles to read.
    # If this number is known in advance it can be given as an input parameter.
    if isnothing(read_positions)

        # this only works if a filter function is provided!
        if isnothing(filter_function)
            error("Please provide either a dictionary with read positions or a filter function!")
        end

        # check how many particles fulfill the filter criterion and find the positions in the blocks.
        # -> needed for array pre-allocation to speed up read-in.
        read_positions = find_read_positions(snap_base, filter_function, verbose=verbose)
    end

    if verbose
        @info "Allocating storage dictionary..."
        flush(stdout)
        flush(stderr)
        t1 = time_ns()
    end

    # select current file
    filename = select_file(snap_base, 0 )

    # read info block
    snap_info = read_info(filename)

    # check if mass block is supposed to be read
    read_mass = !isnothing(findfirst(blocks .== "MASS")) ? true : false

    # check if all blocks are present
    blocks, no_mass_block = check_blocks(filename, blocks, parttype)

    # pre-allocate all data arrays in a dictionary
    d = allocate_data_dict(blocks, read_positions["N_part"], snap_info, no_mass_block)
    
    if verbose
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
        flush(stdout)
        flush(stderr)
    end

    # store the number of particles that have been read
    N_read = 0

    # get the relevant files from the dictionary
    key_entries = collect(keys(read_positions))
    files = key_entries[ key_entries .!= "N_part"]

    # sort the files to always read particles in the same order
    sort!(files)

    # allocate storage arrays for every subfile to read
    # -> needed for multithreading
    num_files = length(files)
    filenames = Vector{String}(undef, num_files)
    headers = Vector{SnapshotHeader}(undef, num_files)
    infos = Vector{Vector{InfoLine}}(undef, num_files)
    block_positions = Vector{Dict{String, Integer}}(undef, num_files)
    indices = Vector{Vector{Int64}}(undef, num_files)
    n_to_reads = Vector{Vector{Int64}}(undef, num_files)
    nreads = Vector{Int64}(undef, num_files)

    # read data for multithreading
    for i ∈ 1:length(files)
        
        # select filename of subfile
        filename = select_file(snap_base, files[i])
        filenames[i] = filename
        # read header block
        headers[i] = read_header(filename)
        # read info block
        infos[i] = read_info(filename)
        # find the location of the blocks in the current file
        block_positions[i] = get_block_positions(filename)
        # store indices of where to start reading
        indices[i] = read_positions[files[i]]["index"]
        # store number of particles to read per index
        n_to_reads[i] = read_positions[files[i]]["n_to_read"]
        # store number of particles that have been read until now
        nreads[i] = N_read
        # count up total particles read until now
        N_read += sum(read_positions[files[i]]["n_to_read"])
    end

    if verbose
        @info "Reading $N_read particles accross $(length(files)) snapshots..."
        flush(stdout)
        flush(stderr)
        t1 = time_ns()
    end

    # to show read progress
    p = Progress(num_files)

    # threaded IO over relevant subfiles
    @threads for (filename, nread, n_to_read, index, block_position, info, h) ∈ collect(zip(filenames, nreads, n_to_reads, indices, block_positions, infos, headers))
       
        # threaded IO over blocks
        @threads for block ∈ blocks

            # open filestream
            f = open(filename, "r")

            block_info = get_requested_info(info, block)

            # store the relevant data in the `d[block]` array
            read_block!(d[block], f, index,
                nread, n_to_read,
                parttype=parttype,
                block_position=block_position[block],
                info=block_info, h=h)

            close(f)
        end

        if verbose
            next!(p)
            flush(stdout)
            flush(stderr)
        end
    end

    if verbose
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
        flush(stdout)
        flush(stderr)
    end

    if read_mass && no_mass_block
        # read header block
        h = read_header(snap_base)
        d["MASS"] = h.massarr[parttype+1] .* ones(Float32, read_positions["N_part"])
    end

    return d
end

"""
    read_blocks_filtered_hdf5( snap_base::String, blocks::Array{String};
                                filter_function::Union{Function, Nothing}=nothing, 
                                read_positions::Union{Dict, Nothing}=nothing, 
                                parttype::Integer=0, verbose::Bool=true )

Reads the specified blocks from all distributed files where particles pass the `filter_function`, or are given by a `Dict` of `read_positions`.
For `read_positions` please see [`find_read_positions`](@ref).
"""
function read_blocks_filtered_hdf5(snap_base::String, blocks::Array{String};
                              filter_function::Union{Function, Nothing}=nothing, 
                              read_positions::Union{Dict, Nothing}=nothing, 
                              parttype::Integer=0, verbose::Bool=true )

    # default behaviour is to find the particles to read.
    # If this number is known in advance it can be given as an input parameter.
    if isnothing(read_positions)

        # this only works if a filter function is provided!
        if isnothing(filter_function)
            error("Please provide either a dictionary with read positions or a filter function!")
        end

        # check how many particles fulfill the filter criterion and find the positions in the blocks.
        # -> needed for array pre-allocation to speed up read-in.
        read_positions = find_read_positions(snap_base, filter_function, verbose=verbose)
    end

    if verbose
        @info "Allocating storage dictionary..."
        flush(stdout)
        flush(stderr)
        t1 = time_ns()
    end

    # select current file
    filename = select_file(snap_base, 0 )

    # check if mass block is supposed to be read
    read_mass = !isnothing(findfirst(blocks .== "Masses")) ? true : false

    # check if all blocks are present
    blocks, no_mass_block = check_blocks(filename, blocks, parttype)

    # pre-allocate all data arrays in a dictionary
    d = allocate_data_dict(blocks, read_positions["N_part"], 
                            h5open(snap_base * ".0.hdf5", "r")["PartType$parttype"], 
                            no_mass_block)
    
    if verbose
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
        flush(stdout)
        flush(stderr)
    end

    # store the number of particles that have been read
    N_read = 0

    # get the relevant files from the dictionary
    key_entries = collect(keys(read_positions))
    files = key_entries[ key_entries .!= "N_part"]

    # sort the files to always read particles in the same order
    sort!(files)


    if verbose
        @info "Reading $(read_positions["N_part"]) particles accross $(length(files)) snapshots..."
        flush(stdout)
        flush(stderr)
        t1 = time_ns()
    end

    # to show read progress
    p = Progress(length(files))
        
    # number of particles read so far
    n_read = 0

    for file ∈ files
        
        # construct indices at which the file should be read
        indices = [read_positions[file]["index"][i]+1:read_positions[file]["index"][i]+read_positions[file]["n_to_read"][i] 
                    for i ∈ 1:length(read_positions[file]["index"])]
        indices = vcat(collect.(indices)...)

        # nuumber of particles to read in this file
        n_to_read = length(indices)
        
        # open file
        f = h5open(select_file(snap_base, file), "r")
        # threaded IO over blocks
        @threads for block ∈ blocks
            if ndims(f["PartType$parttype"][block]) == 1
                d[block][1+n_read:n_read+n_to_read] = HDF5.readmmap(f["PartType$parttype"][block])[indices]
            else
                ndim = size(f["PartType$parttype"][block])[1]
                d[block][1:ndim,1+n_read:n_read+n_to_read] = HDF5.readmmap(f["PartType$parttype"][block])[1:ndim, indices]
            end
        end

        n_read += n_to_read

        # close file
        close(f)
        if verbose
            next!(p)
            flush(stdout)
            flush(stderr)
        end

    end

    if verbose
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
        flush(stdout)
        flush(stderr)
    end

    if read_mass && no_mass_block
        # read header block
        h = read_header(snap_base)
        d["MASS"] = h.massarr[parttype+1] .* ones(Float32, read_positions["N_part"])
    end

    return d
end