"""
    check_block_format(filename)

Checks which binary format the file is in.
"""
function check_snapshot_format(filename)

    # ToDo: HDF5!

    f = open(filename)
    blocksize = read(f, Int32)
    close(f)

    if blocksize == 8
        # snapshot format 2
        return 2
    elseif blocksize == 256
        # snapshot format 1
        return 1
    else
        error("File is not in Gadget binary format!")
    end
end

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

        # check for up-to 10x integer overflow - everything beyond that is ridiculous anyway
        for N_check = 1:10
            # sum up actual block size
            blocksize_before_fix = blocksize_before + N_check*4294967296
            seek(f,position_before+blocksize_before_fix+12)
            # read potential blocksize at the end of the block
            blocksize_after_fix = read(f,UInt32) + N_check*4294967296
            
            # check if blocksizes match
            if blocksize_before_fix == blocksize_after_fix
                return blocksize_before_fix
            end
        end

        error("There is an issue with the snapshot:\n
                   Blocksize_before = $blocksize_before_fix\n
                   Blocksize_after  = $blocksize_after_fix")
    end
end

"""
    check_blocksize_format1(f::IOStream, position_before::Integer, blocksize_before::Integer)

Checks for integer overflow in the size of the block.
"""
function check_blocksize_format1(f::IOStream, position_before::Integer, blocksize_before::Integer)

    seek(f,position_before+blocksize_before)
    blocksize_after = read(f,UInt32)

    # if the numbers are the same everything works as it should
    if blocksize_before == blocksize_after
        return blocksize_before
    else # compensate for integer overflow

        # check for up-to 10x integer overflow - everything beyond that is ridiculous anyway
        for N_check = 1:10
            # sum up actual block size
            blocksize_before_fix = blocksize_before + N_check*4294967296
            seek(f,position_before+blocksize_before_fix)
            # read potential blocksize at the end of the block
            blocksize_after_fix = read(f,UInt32) + N_check*4294967296
            
            # check if blocksizes match
            if blocksize_before_fix == blocksize_after_fix
                return blocksize_before_fix
            end
        end

        error("There is an issue with the snapshot:\n
                   Blocksize_before = $blocksize_before_fix\n
                   Blocksize_after  = $blocksize_after_fix")
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

    # if hdf5 file use seperate function
    if HDF5.ishdf5(filename)
        return check_blocks_hdf5(filename, blocks, parttype)
    end

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
    check_blocks(snap_base::String, blocks::Array{String})

Check if all requested blocks are present.
"""
function check_blocks_hdf5(filename::String, blocks::Array{String}, parttype::Integer)

    # check if requested blocks are present
    no_mass_block = false

    f = h5open(filename, "r")
    blocks_in_file = keys(f["PartType$parttype"])

    for blockname ∈ blocks
        if !block_present(filename, blockname, blocks_in_file)
            if blockname == "Masses"
                no_mass_block = true
                deleteat!(blocks, findfirst(blocks .== "Masses"))
            else
                error("Block $blockname not present!")
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

    # check if file is present
    if isfile(filebase)
        return filebase
    else
        # check if filebase is a hdf5 file
        if isfile(filebase * ".hdf5")
            return filebase * ".hdf5"
        else
            # check if subfile is present
            filename = filebase * ".$filenum"
            if isfile(filename)
                return filename
            else
                # check if subfile is hdf5 file
                if isfile(filename * ".hdf5")
                    return filename * ".hdf5"
                else
                    error("File $filename not present!")
                end
            end # subfile            
        end # filebase hdf5
    end # filebase correct file
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
