"""
    get_block_positions(filename::String; 
                        snap_format::Integer=2)

Returns a dictionary with the starting positions of all blocks in a snapshot in bits.

# Example
```julia
# Format 2 -> default
block_positions = get_block_positions(filename)
block_positions["POS"]
# Format 1
block_positions = get_block_positions(filename)
block_positions[1]
```
"""
function get_block_positions(filename::String)

    # get snapshot format
    snap_format = check_snapshot_format(filename)

    if snap_format == 1
        return get_block_positions_format1(filename)
    elseif snap_format == 2
        return get_block_positions_format2(filename)
    end
end
"""
    get_block_positions_format2(filename::String)

Returns a dictionary with the starting positions of all blocks in a snapshot in bits.

# Example 
```julia
block_positions = get_block_positions(filename)
block_positions["POS"]
```
"""
function get_block_positions_format2(filename::String)

    f = open(filename)
    blocksize = read(f, Int32)

    # only works for snap format 2
    if blocksize != 8
        error("Block search not possible - use snap_format 2!")
    end

    # allocate data dict
    d = Dict{String, Integer}()

    while eof(f) != true

        # read block name
        blockname = read_bockname(f)

        p = position(f)

        seek(f,p+8)

        skipsize = read(f, UInt32)

        skipsize = check_blocksize(f, p, skipsize)

        # store blockname and position in Dict
        d[blockname] = p+12

        # skip to the next name block
        seek(f,p+skipsize+20)

    end

    close(f)

    return d
end

"""
    get_block_positions_format1(filename::String)

Returns a dictionary with the starting positions of all blocks in a snapshot in bits.

# Example 
```julia
block_positions = get_block_positions(filename)
block_positions["POS"]
```
"""
function get_block_positions_format1(filename::String)

    f = open(filename)

    # read the size of the data block -> should always be 264!
    blocksize = read(f, Int32)

    if blocksize != 256
        error("Reading a block by number only works for Format 1!")
    end

    # allocate data dict
    d = Dict{Integer, Integer}()

    # skip to the end of the header
    seek(f, 264)

    # define block number
    i_block = 1

    while eof(f) != true

        # read the size of the data block -> should always be 264!
        blocksize = read(f, UInt32)

        # beginning of data in block
        p = position(f)

        # check for integer overflow in large blocks
        blocksize = check_blocksize_format1(f, p, blocksize)

        # store blockname and position in Dict
        d[i_block] = p

        # skip to the next name block
        seek(f,p+blocksize+4)

        # count up block number
        i_block += 1
    end

    close(f)

    return d
end


"""
    check_block_position(filename::String, blockname::String)

Helper function to find read positions.
"""
function check_block_position(filename::String, blockname::String)
    # read positions of all blocks
    block_positions_dict = get_block_positions(filename)

    # check if the requested block is present
    if haskey(block_positions_dict, blockname)
        return block_positions_dict[blockname], blockname == "MASS"
    else
        if blockname != "MASS"
            # if the block is not present we need error handling!
            error("Requested block $blockname not present!")
        else
            return -1, true
        end
    end
end


"""
    check_block_position(filename::String, blocknum::Integer)

Helper function to find read positions.
"""
function check_block_position(filename::String, blocknum::Integer)
    # read positions of all blocks
    block_positions_dict = get_block_positions_format1(filename)

    # check if the requested block is present
    if haskey(block_positions_dict, blocknum)
        return block_positions_dict[blocknum], blocknum == 4
    else
        # mass block is always number 4
        if blocknum != 4
            # if the block is not present we need error handling!
            error("Requested block $blocknum not present!")
        else
            return -1, true
        end
    end
end