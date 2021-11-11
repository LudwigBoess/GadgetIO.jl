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
function read_particles_in_box_peano(filename::String, blocks::Vector{String},
                               corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                               parttype::Integer=0, verbose::Bool=true)


    if verbose
        println()
        @info "Running on $(nthreads()) threads"
    end

    filebase = filename
    file_key_index = filebase * ".key.index"

    # if the snapshot does not exist it may be split into multiple files
    if !isfile(filebase)

        if verbose
            @info "File: $filebase not found, looking for sub-files."
        end

        # try reading the first of the distributed snapshots
        filename = select_file(filebase, 0)

        h = head_to_obj(filename)

        if verbose
            @info "$(h.num_files) sub-files found."
        end

        nfiles = h.num_files
    else
        nfiles = 1
        h = head_to_obj(filename)
    end

    blocks, no_mass_block = check_blocks(filename, blocks)

    # read info blocks once here
    snap_info = read_info(filename)
    key_info  = read_info(filename * ".key")

    if verbose
        @info "All requested blocks present!"
        @info "Checking for .key files..."
    end

    # check if key files are present
    file_key = filename * ".key"
    if !block_present(file_key, "KEY")
        error("No .key file present!")
    end

    if verbose
        @info ".key files found!"
        @info "Calculating peano-hilbert keys..."
        t1 = Dates.now()
    end

    # first read the header
    h_key = read_keyheader(file_key)

    # get a list of the required peano-hilbert keys
    keylist = get_keylist(h_key, corner_lowerleft, corner_upperright)

    if verbose
        t2 = Dates.now()
        @info "$(size(keylist,1)) Peano-Hilbert keys found. Took: $(t2 - t1)"
        @info "Looking for relevant files..."
        t1 = Dates.now()
    end

    # find relevant files
    files = find_files_for_keys(filebase, nfiles, keylist)
    
    N_files = size(files,1)
    if verbose
        t2 = Dates.now()
        @info "$N_files files found. Took: $(t2 - t1)"

        @info "Searching read positions..."
        println()
        t1 = Dates.now()
    end

    # find all the positions where to read data
    file_offset_key, file_part_per_key, file_block_positions = find_read_positions( files, filebase, blocks, 
                                                                                    parttype, keylist, key_info, 
                                                                                    verbose)

    N_to_read = 0

    @inbounds for i = 1:N_files
        N_to_read += sum(file_part_per_key[i])
    end

    if verbose
        t2 = Dates.now()
        println()
        @info "Positions read. Took: $(t2 - t1)"
        println()
        @info "Reading $N_to_read particles..."
    end

    # prepare dictionary for particle storage
    d = allocate_data_dict(blocks, N_to_read, snap_info, no_mass_block)

    if verbose
        @info "Reading Blocks..."
        t1 = Dates.now()
    end

    n_read = 0

    for i = 1:N_files

        # select current file
        filename = select_file(filebase, files[i])

        # read header
        h = read_header(filename)

        # no particles of parttype in the file
        if h.npart[parttype+1] == 0
            continue
        end

        # read info block
        snap_info = read_info(filename)


        # read blocks in parallel
        @threads for j = 1:size(blocks,1)

            block_info = snap_info[getfield.(snap_info, :block_name) .== blocks[j]][1]

            # add offset of particle types that should not be read
            offset = 0
            for k=1:size(h.npart,1)
                if (block_info.is_present[k] > 0) & (h.npart[k] > 0) & ( k < parttype + 1)
                    offset += h.npart[k]
                end
            end

            # reads data into the dictionary and counts up n_read
            read_block_with_offset!(d[blocks[j]], n_read, filename, 
                                    file_block_positions[i][blocks[j]],
                                    block_info, offset, file_offset_key[i],
                                    file_part_per_key[i])

        end # loop over blocks

        
        n_read += sum(file_part_per_key[i])

        @info "Read $n_read / $N_to_read particles"

    end # for i = 1:size(files,1)

    # finally construct masses of no mass block present
    if no_mass_block
        d["MASS"] = h.massarr[parttype+1] .* ones(Float32, N_to_read)
    end

    if verbose
        t2 = Dates.now()
        @info "Blocks read. Took: $(t2 - t1)"
    end

    return d
end


"""
    read_particles_in_box(filename::String, blocks::String,
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
        d = read_blocks_over_all_files(filename, blocks, filter_function = filter_function, parttype = parttype, verbose = verbose )
    end

    return d
end

"""
    read_particles_in_box(filename::String, blocks::String,
                          corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                          parttype::Integer=0, verbose::Bool=true,
                          use_keys::Bool=true)

Like `read_particles_in_box` but for a single block. Returns the block as an array.

See also: [`read_particles_in_box`](@ref)
"""
function read_particles_in_box(filename::String, blocks::String,
                               corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real};
                               parttype::Integer=0, verbose::Bool=true,
                               use_keys::Bool=true)

    d = read_particles_in_box(filename, [blocks], corner_lowerleft, corner_upperright, 
                              parttype=parttype, verbose=verbose, use_keys=use_keys)

    return d[blocks]
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
    read_particles_in_box(filename::String, blocks::String,
                          center_pos, radius;
                          parttype::Integer=0, verbose::Bool=true)

Like `read_particles_in_volume` but for a single block. Returns the block as an array.

See also: [`read_particles_in_volume`](@ref)
"""
function read_particles_in_volume(filename::String, blocks::String,
                                  center_pos, radius;
                                  parttype::Integer=0, verbose::Bool=true,
                                  use_keys::Bool=true)

    d = read_particles_in_volume(filename, [blocks], center_pos, radius,
                                 parttype=parttype, verbose=verbose,
                                 use_keys=use_keys)

    return d[blocks]
end

