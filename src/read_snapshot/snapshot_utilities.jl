
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

        name = Char.(read!(f, Array{Int8,1}(undef,4)))
        blockname = String(name)

        blockname = strip(blockname)

        push!(blocks, blockname)

        p = position(f)
        seek(f,p+8)

        skipsize = read(f, Int32)

        p = position(f)
        seek(f,p+skipsize+8)

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
    read_info(filename; verbose::Bool=false)

Reads the info block of a snapshot and returns the information in an array of `Info_Line` types.
If `verbose=true` the blocknames are also printed to console.
"""
function read_info(filename::String; verbose::Bool=false)

    f = open(filename)
    seek(f,4)

    while eof(f) != true

        name = Char.(read!(f, Array{Int8,1}(undef,4)))
        blockname = String(name)

        p = position(f)
        seek(f,p+8)

        skipsize = read(f, Int32)

        if blockname == "INFO"
            n_blocks = Int(skipsize/40) # one info line is 40 bytes
            arr_info = Array{Info_Line,1}(undef,n_blocks)

            for i = 1:n_blocks
                arr_info[i] = read_info_line(f)
            end # for

            close(f)

            if verbose == true
                println("Found Info block.\nEntries are:\n")
                for i ∈ 1:length(arr_info)
                    println(i, " - ", arr_info[i].block_name)
                end # for
            end # verbose

            return arr_info

        else
            p = position(f)
            seek(f,p+skipsize+8)
        end # if blockname == "INFO"

    end # while eof(f) != true

    println("No info block present!")

    return 1
end

function read_info_line(f)

    # block name consists of 4 C Chars.
    name = Char.(read!(f, Array{Int8,1}(undef,4)))
    blockname = String(name)

    block_name = strip(blockname)

    # the datatype is stored as a 8 char word.
    letters = read!(f, Array{Int8,1}(undef,8))
    letters = Char.(letters)
    letters = string.(letters)

    data_type = ""
    for i = 1:length(letters)
        data_type *= letters[i]
    end

    # erase whitespace and make lower case for comparison
    data_type = lowercase(strip(data_type))

    # assign data types
    if data_type == "float" ||
       data_type == "floatn"
            dt = Float32
    elseif data_type == "long"
            dt = UInt32
    elseif data_type == "llong"
            dt = UInt64
    elseif data_type == "double" ||
           data_type == "doublen"
            dt = Float64
    end

    # array dimensions are stored as Int32
    n_dim = read(f,Int32)

    # read the information for which particles the block is relevant.
    # e.g. for gas particles: [ 1, 0, 0, 0, 0, 0 ]
    is_present = read!(f, Array{Int32,1}(undef,6))

    # construct the info line struct and return it.
    return Info_Line(block_name, dt, n_dim, is_present)
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
function check_blocks(filename::String, blocks::Array{String})

    # check if requested blocks are present
    no_mass_block = false
    blocks_in_file = print_blocks(filename, verbose=false)
    for blockname in blocks
        if !block_present(filename, blockname, blocks_in_file)
            if blockname == "MASS"

                no_mass_block = true
                deleteat!(blocks, "MASS")
            else
                error("Block $blockname not present!")
            end
        end
    end

    return blocks, no_mass_block

end


"""
    get_block_positions(filename::String)

Returns a dictionary with the starting positions of all blocks in a snapshot in bits.
"""
function get_block_positions(filename::String)

    f = open(filename)
    blocksize = read(f, Int32)

    # only works for snap format 2
    if blocksize != 8
        error("Block search not possible - use snap_format 2!")
    end

    p = position(f)
    seek(f,4)

    blocks = Vector{String}(undef, 0)
    pos    = Vector{Int64}(undef, 0)


    while eof(f) != true

        # read block name
        name = Char.(read!(f, Array{Int8,1}(undef,4)))
        blockname = String(name)

        blockname = strip(blockname)

        # store block name in array
        push!(blocks, blockname)

        read(f, Int32)
        read(f, Int32)

        skipsize = read(f, Int32)

        p = position(f)

        push!(pos, p)

        seek(f,p+skipsize+8)

    end

    close(f)

    d = Dict( blocks[i] => pos[i] for i = 1:length(blocks))

    return d
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
    allocate_data_dict( blocks::Array{String}, N_to_read::Integer, 
                        snap_info::Array{Info_Line}, no_mass_block::Bool )

Helper function to allocate the data `Dict`.
"""
function allocate_data_dict(blocks::Array{String}, N_to_read::Integer, 
                            snap_info::Array{Info_Line}, no_mass_block::Bool)

    # prepare dictionary for particle storage, this looks super ugly...
    d = Dict(blocks[i] => Array{snap_info[getfield.(snap_info, :block_name) .== blocks[i]][1].data_type,2}(undef, N_to_read,
            snap_info[getfield.(snap_info, :block_name) .== blocks[i]][1].n_dim) for i = 1:length(blocks))

    # allocate mass array, if it's not in a block
    if no_mass_block
        d["MASS"] = Array{Float32,2}(undef, N_to_read, 1)
    end

    return d
end
