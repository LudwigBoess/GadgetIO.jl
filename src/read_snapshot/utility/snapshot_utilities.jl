"""
    check_blocksize(f::IOStream, position_before::Integer, blocksize_before::Integer)

Checks for integer overflow in the size of the block.
"""
function check_blocksize(f::IOStream, position_before::Integer, blocksize_before::Integer)

    seek(f,position_before+blocksize_before+12)
    blocksize_after = read(f,UInt32)

    # if the numbers are the same everything works as it should
    if blocksize_before == blocksize_after
        return blocksize_before
    else # compensate for integer overflow
        blocksize_before_fix = blocksize_before + 4294967296
        seek(f,position_before+blocksize_before_fix+12) 
        blocksize_after_fix = read(f,UInt32) + 4294967296
        
        if blocksize_before_fix == blocksize_after_fix
            return blocksize_before_fix
        else
            error("There is an issue with the snapshot:\n
                   Blocksize_before = $blocksize_before_fix\n
                   Blocksize_after  = $blocksize_after_fix")
        end
    end
end

"""
    read_bockname(f::IOStream)

Reads the name of the Format 2 block.
"""
function read_bockname(f::IOStream)
    name = Char.(read!(f, Array{Int8,1}(undef,4)))
    return strip(String(name))
end

"""
    print_blocks(filename::String; verbose::Bool=true)

Reads the block names of blocks in a snapshot and returns them in an array.
Outputs them to console if `verbose=true`
"""
function print_blocks(filename::String; verbose::Bool=true)

    f = open(filename)
    blocksize = read(f, Int32)

    if blocksize != 8
        error("Block search not possible - use snap_format 2!")
    end

    p = position(f)
    seek(f,4)

    blocks = []

    while eof(f) != true

        # read block name
        blockname = read_bockname(f)

        # store block name in array
        push!(blocks, blockname)

        p = position(f)

        seek(f,p+8)

        skipsize = read(f, UInt32)

        skipsize = check_blocksize(f, p, skipsize)

        seek(f,p+skipsize+20)

    end

    if verbose
        println("Found blocks: ")
        for block ∈ blocks
            println(block)
        end
    end

    close(f)

    return String.(blocks)
end


"""
    block_present(filename::String, blockname::String, blocks::Vector{String}=[""])

Checks if a given block is present in the snapshot file, or in the supplied `blocks` Vector.
"""
function block_present(filename::String, blockname::String, blocks::Vector{String}=[""])

    # if no blocks supplied, read in block information
    if blocks == [""]
        blocks = print_blocks(filename, verbose=false)
    end

    # loop over all blocks
    for block ∈ blocks
        if block == blockname
            return true
        end
    end

    return false
end


"""
    check_blocks(snap_base::String, blocks::Array{String})

Check if all requested blocks are present.
"""
function check_blocks(filename::String, blocks::Array{String}, parttype::Integer)

    # check if requested blocks are present
    no_mass_block = false
    blocks_in_file = print_blocks(filename, verbose=false)
    # read info block
    snap_info = read_info(filename)

    for blockname ∈ blocks
        if !block_present(filename, blockname, blocks_in_file)
            if blockname == "MASS"
                no_mass_block = true
                deleteat!(blocks, findfirst(blocks .== "MASS"))
            else
                error("Block $blockname not present!")
            end
        else
            if blockname == "MASS" && iszero(snap_info[getfield.(snap_info, :block_name) .== "MASS"][1].is_present[parttype+1])
                no_mass_block = true
                deleteat!(blocks, findfirst(blocks .== "MASS"))
            end
        end
    end

    return blocks, no_mass_block

end



"""
    select_file(filebase::String, filenum::Integer)

Checks if the requested files exists and returns the path.
"""
function select_file(filebase::String, filenum::Integer)

    if !isfile(filebase)
        filename = filebase * ".$filenum"
        if !isfile(filename)
            error("File $filename not present!")
        else
            return filename
        end
    else
        return filebase
    end
end

""" 
    shift_across_box_border(x::Real, x_halo::Real, boxsize::Real, boxsize_half::Real)

Shift coordinate `x` across the box border if the zero coordinate `x₀` is on the other side.
"""
@inline function shift_across_box_border(x::Real, x₀::Real, boxsize::Real, boxsize_half::Real)
    if x - x₀ > boxsize_half
        return x - boxsize
    elseif x₀ - x > boxsize_half
        return x + boxsize
    end
    return x
end
