"""
    check_block_position(filename::String, blockname::String)

Helper function to find read positions.
"""
function check_block_position(filename::String, blockname::String)
    # read positions of all blocks
    block_positions_dict = get_block_positions(filename)

    # check if the requested block is present
    if haskey(block_positions_dict, blockname)
        return block_positions_dict[blockname], false
    else
        if blockname != "MASS"
            # if the block is not present we need error handling!
            error("Requested block $blockname not present!")
        else
            return -1, true
        end

    end
end

