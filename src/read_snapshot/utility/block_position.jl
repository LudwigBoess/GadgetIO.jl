"""
    get_block_positions(filename::String)

Returns a dictionary with the starting positions of all blocks in a snapshot in bits.

# Example 
```julia
block_positions = get_block_positions(filename)
block_positions["POS"]
```
"""
function get_block_positions(filename::String)

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

