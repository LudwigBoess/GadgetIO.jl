"""
    read_block(filename::String, block_id::Union{String, Integer};
                parttype::Integer=0,
                block_position::Integer=-1,
                info::Union{Nothing,InfoLine}=nothing,
                h::Union{Nothing,SnapshotHeader}=nothing,
                offset=0, n_to_read=-1)

Reads a block in a snapshot with `block_id` as a given name or position in the file. Defaults to reading gas particles.
Block Names are case sensitive.

# Keyword Arguments
- `parttype`: Which particle type to read (0-indexed)
    - `0`: Gas (default)
    - `1`: Dark Matter
    - `2`: Boundary
    - `3`: Bulge
    - `4`: Stars 
    - `5`: Black Holes
    - `-1`: All particle types
- `block_position`: Position of the block in the data file `[bytes]`. Can be used to speed up IO.
- `info`: `InfoLine` for the given block. Can be used to speed up IO, or if no `INFO` block is present.
- `h`: `SnapshotHeader` for given file. Can be used to speed up IO.
- `offset`: Adds an offset to start reading the block at a later point. Can be used for sub-IO.
- `n_to_read`: How many particles to read. Can be used for sub-IO.

# Examples (Format 2)
In case an `INFO` block is present:
```julia
gas_pos = read_block(filename, "POS", parttype=0)
```
If not you need to supply your own `InfoLine`
```julia
pos_info = InfoLine("POS", Float32, 1, [1, 1, 1, 1, 1, 1])
gas_pos = read_block(filename, "POS", info=pos_info, parttype=0)

# Examples (Format 1)
In case you want to read the default blocks:
```julia
gas_pos = read_block(filename, 1, parttype=0)
```
If not you need to supply your own `InfoLine`
```julia
pos_info = InfoLine("", Float32, 1, [1, 1, 1, 1, 1, 1])
gas_pos = read_block(filename, 1, info=pos_info, parttype=0)
```
"""
function read_block(filename::String, block_id::Union{String, Integer};
                    parttype::Integer=0,
                    block_position::Integer=-1,
                    info::Union{Nothing,InfoLine}=nothing,
                    h::Union{Nothing,SnapshotHeader}=nothing,
                    offset::Integer=0, n_to_read::Integer=-1)


    # handle reading all particle types seperately to avoid clutter in this function
    if parttype == -1
        return read_all_parttypes(filename, block_id,
            block_position=block_position, info=info, h=h)
    end

    if isnothing(h)
        # read header - super fast and needed for flexibility
        h = read_header(filename)
    end

    # define default behaviour here
    n_read_io = true

    # if number of entries to read is not specified
    if n_to_read == -1
        n_read_io = false

        if isfile(filename)
            # if there are no subfiles or we are already reading a subfile
            # -> read local number of particles

            # if specific particle species needs to be read
            n_to_read = h.npart[parttype+1]
            # we have calculated the number of particles to read now
            n_read_io = true
        else
            # if there are sub-files
            # -> read all entries

            # calculate total number of particles to read
            n_to_read = get_total_particles(h, parttype)
        end
    end

    # check if particle type is present
    if iszero(n_to_read)
        error("Particle Type $parttype not present in simulation!")
    end

    if isfile(filename)
        num_files = 1
    else
        num_files = h.num_files
    end

    # check if info is present
    was_info_given = !isnothing(info)
    if isnothing(info)
        info = check_info(filename, block_id)
    end

    # store if block position has been provided
    block_position_io = false

    # check if block position has been provided
    if block_position != -1
        block_position_io = true
    end

    # allocate storage here
    block = allocate_data_array(info, n_to_read)

    # we have not read any data yet
    nread = 0
    offset_rest = offset # remaining offset (relevant when the offset is larger than a file is long)

    # loop over all sub-files
    # -> irrelevant if already a subfile
    filenames = Vector{String}(undef, num_files)
    headers = Vector{SnapshotHeader}(undef, num_files)
    infos = Vector{InfoLine}(undef, num_files)
    nreads = Vector{Int64}(undef, num_files)
    n_to_reads = Vector{Int64}(undef, num_files)
    block_positions = Vector{Int64}(undef, num_files)
    offs = Vector{Int64}(undef, num_files)
    for file ∈ 0:num_files-1
        
        # read local filename
        _filename = select_file(filename, file)

        if n_to_read > nread
            h = read_header(_filename)

            # get number of particles to read from local file, if not given
            n_to_read_file = h.npart[parttype+1]
            if n_to_read_file > offset_rest
                n_to_read_file -= offset_rest
                off = offset_rest
                offset_rest = 0
            else
                off = 0
                offset_rest -= n_to_read_file
                n_to_read_file = 0
            end
            n_to_read_file = min(n_to_read_file, n_to_read - nread)
        else
            # do not read anything anymore
            n_to_read_file = 0
        end

        if iszero(n_to_read_file)
            ind = file + 1
            filenames[ind] = _filename
            headers[ind] = h
            infos[ind] = InfoLine(block_id, eltype(block), 0, zeros(Int32, 6))
            nreads[ind] = nread
            n_to_reads[ind] = n_to_read_file
            block_positions[ind] = 0
            offs[ind] = 0
            continue
        end

        # read info individually for the file if no info block was passed
        if !was_info_given
            info = check_info(_filename, block_id)
        end

        # find block position, if not given
        if !block_position_io
            # read block position, if not given
            block_position, mass_block = check_block_position(_filename, block_id)

            # check if the block is even present for the requested particle
            if !iszero(h.npart[parttype+1]) && (iszero(info.is_present[parttype+1]) && !mass_block)
                error("Requested block $block_id not present for particle type $(parttype)!")
            end

            # if the mass block is not present for the current particle type it's stored in the header 
            if !iszero(h.npart[parttype+1]) && mass_block && !iszero(h.massarr[parttype+1])
                # assign mass from header
                block .= h.massarr[parttype+1]
                # we can exit the function now
                return block
            end
        end

        # store info for current subfile
        ind = file + 1
        filenames[ind] = _filename
        headers[ind] = h
        infos[ind] = info
        nreads[ind] = nread
        n_to_reads[ind] = n_to_read_file
        block_positions[ind] = block_position
        offs[ind] = off

        # count up number of particles read including this subfile
        nread += n_to_read_file
    end

    # threaded IO over all subfiles
    Threads.@threads for (_filename, nread, n_to_read, block_position, off, info, h) ∈ collect(zip(filenames, nreads, n_to_reads, block_positions, offs, infos, headers))
        if iszero(n_to_read)
            continue
        end

        f = open(_filename, "r")
        read_block!(block, f, off, nread, n_to_read;
                    parttype, block_position, info, h)

        close(f)
    end

    return block
end


"""
    read_block!(a::AbstractArray, f::IOStream,
                offset::Integer,
                nread::Integer, n_to_read::Integer;
                parttype::Integer,
                block_position::Integer,
                info::InfoLine,
                h::SnapshotHeader)

Read part of a block from a given file stream into a pre-allocated array.
"""
function read_block!(a::AbstractArray, f::IOStream,
                    offset::Integer=0,
                    nread::Integer=0, n_to_read::Integer=-1;
                    parttype::Integer,
                    block_position::Integer,
                    info::InfoLine,
                    h::SnapshotHeader)

    # number of bits in data_type
    len = sizeof(info.data_type) * info.n_dim

    # offset due to particles before requested one
    part_offset = 0

    # count up offset for relevant blocks
    if parttype > 0
        for i = 1:parttype
            part_offset += h.npart[i] * info.is_present[i]
        end
    end


    # jump to position of particle type in relevant block
    seek(f, block_position+(part_offset+offset)*len)

    # note the Int(...) are necessary to enable the reading into the SubArray returned by @view
    if info.n_dim == 1
        read!(f, @view(a[Int(nread+1):Int(nread+n_to_read)]))
    else
        read!(f, @view(a[:, Int(nread+1):Int(nread+n_to_read)]))
    end

    return a

end

"""
    read_block!(a::AbstractArray, f::IOStream,
                offset::Array{<:Integer},
                nread::Integer, n_to_read::Array{<:Integer};
                parttype::Integer,
                block_position::Integer,
                info::InfoLine,
                h::SnapshotHeader)

Read part of a block from a given file stream into a pre-allocated array.
"""
function read_block!(a::AbstractArray, f::IOStream,
                    offset::Array{<:Integer},
                    nread::Integer, n_to_read::Array{<:Integer};
                    parttype::Integer,
                    block_position::Integer,
                    info::InfoLine,
                    h::SnapshotHeader)

    # store local variables
    nread_local  = nread

    # loop over all offsets
    for i = 1:length(offset)

        read_block!(a, f, offset[i], nread_local, n_to_read[i];
                    parttype, block_position, info, h)

        # store locally read particles
        nread_local  += n_to_read[i]
    end

    a
end


"""
    read_all_parttypes( filename::String, block_id::String;
                        block_position::Integer=-1,
                        info::Union{Nothing,InfoLine}=nothing,
                        h::Union{Nothing,SnapshotHeader}=nothing)

Reads the requested block for all particle types.
"""
function read_all_parttypes(filename::String, block_id::Union{String, Integer};
                            block_position::Integer=-1,
                            info::Union{Nothing,InfoLine}=nothing,
                            h::Union{Nothing,SnapshotHeader}=nothing)

    # read the header if it's not provided
    if isnothing(h)
        h = read_header(filename)
    end

    # read the info line of not provided
    if isnothing(info)
        info = check_info(filename, block_id)
    end

    # sum up the total number of particles in the simulation
    n_to_read = Vector{Int64}(undef, 6)
    for parttype = 0:5
        if isfile(filename)
            n_to_read[parttype+1] = h.npart[parttype+1]
        else
            # calculate total number of particles to read
            n_to_read[parttype+1] = get_total_particles(h, parttype)
        end
    end

    # allocate storage here
    block = allocate_data_array(info, sum(n_to_read))

    # store number of particles that have been read so far
    nread = 0
    for parttype = 0:5
        # only read the block if particles are actually present
        if !iszero(n_to_read[parttype+1])
            if info.n_dim == 1
                    block[Int(nread + 1):Int(nread + n_to_read[parttype+1])] = read_block(filename, block_id, parttype=parttype,
                                                                    block_position=block_position, info=info, h=h)
            else
                    block[:, Int(nread + 1):Int(nread + n_to_read[parttype+1])] = read_block(filename, block_id, parttype=parttype,
                                                                    block_position=block_position, info=info, h=h)
            end
            nread += n_to_read[parttype+1]
        end
    end

    return block
end
