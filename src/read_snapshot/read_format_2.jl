"""
read_block(filename::String, blockname::String;
           parttype::Integer=0,
           block_position::Integer=-1,
           info::Union{Nothing,InfoLine}=nothing,
           h::Union{Nothing,SnapshotHeader}=nothing,
           offset=0, n_to_read=-1)

Reads a block in a snapshot with given name. Defaults to reading gas particles.
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

# Examples
In case an `INFO` block is present:
```julia
gas_pos = read_block(filename, "POS", parttype=0)
```
If not you need to supply your own `InfoLine`
```julia
pos_info = InfoLine("POS", Float32, 1, [1, 1, 1, 1, 1, 1])
gas_pos = read_block(filename, "POS", info=pos_info, parttype=0)
```
"""
function read_block(filename::String, blockname::String;
                    parttype::Integer=0,
                    block_position::Integer=-1,
                    info::Union{Nothing,InfoLine}=nothing,
                    h::Union{Nothing,SnapshotHeader}=nothing,
                    offset::Integer=0, n_to_read::Integer=-1)

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
            if parttype != -1
                n_to_read = h.npart[parttype+1]
            else
                n_to_read = sum(h.npart)
            end
        else
            # if there are sub-files
            # -> read all entries

            # if specific particle species needs to be read
            if parttype != -1
                n_to_read = get_total_particles(h, parttype)
            else
                n_to_read = 0
                # sum up all particles
                for ptype = 0:5
                    n_to_read += get_total_particles(h, ptype)
                end
            end

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
    if isnothing(info)
        info = check_info(filename, blockname)
    end

    block_position_io = false

    if block_position != -1
        block_position_io = true
    end

    # allocate storage here
    block = allocate_data_array(info, n_to_read)

    # we have not read any data yet
    nread = 0

    # loop over all sub-files
    # -> irrelevant if already a subfile
    for file ∈ 0:num_files-1
        
        # read local filename
        _filename = select_file(filename, file)

        # get number of particles to read from local file, if not given
        if (h.num_files > 1) && (!n_read_io)
            h_internal = read_header(_filename)
            n_to_read  = h_internal.npart[parttype+1]
        end

        # find block position, if not given
        if !block_position_io
            # read block position, if not given
            block_position, mass_block = check_block_position(_filename, blockname)

            if mass_block
                block .= h.massarr[parttype+1]
                return block
            end
        end
        
        f = open(_filename, "r")
        read_block!(block, f, offset, nread, n_to_read;
                    parttype, block_position, info, h)

        close(f)
        nread += n_to_read
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

