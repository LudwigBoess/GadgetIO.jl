"""
    read_block(filename::String, blockname::String;
                                info::InfoLine=InfoLine(),
                                parttype::Integer=-1)

Reads a block in a snapshot with given name. Names are case sensitive.

# Examples
```jldoctest
julia> pos_info = InfoLine("POS", Float32, 1, [1, 1, 1, 1, 1, 1])
[...]
julia> gas_pos = read_block(filename, "POS", info=pos_info, parttype=0)
[...]
```
"""
function read_block(filename::String, blockname::String;
                    parttype::Integer=-1,
                    block_position::Integer=-1,
                    info::Union{Nothing,InfoLine}=nothing,
                    h::Union{Nothing,SnapshotHeader}=nothing)

    # we need to know which particle to read
    if parttype == -1
        error("Please specify particle type!")
    end

    if isnothing(h)
        # read header - super fast and needed for flexibility
        h = head_to_obj(filename)
    end

    # check if particle type is present
    if h.nall[parttype+1] == UInt32(0)
        error("Particle Type $parttype not present in simulation!")
    end

    # check if info is present
    if isnothing(info)
        info = check_info(filename, blockname)
    end

    # check if a block position is supplied
    if block_position == -1

        block_position = check_block_position(filename, blockname)

        if block_position == -1
            # if no mass block is present we can read it from the header
            if blockname == "MASS"
                block = Array{info.data_type,2}(undef,(h.npart[parttype+1], info.n_dim))
                block .= h.massarr[parttype+1]
                return block
            end
        end
    end


    f = open(filename)

    seek(f,block_position)

    for i ∈ 1:size(info.is_present)[1]
        p = position(f)

        if info.is_present[i] == Int32(1)
            if i == (parttype+1)
                block = read_block_data(f, info.data_type, info.n_dim, h.npart[i])
                close(f)
                return block
            else
                seek(f, p + ( sizeof(info.data_type)*info.n_dim*h.npart[i] ))
            end # if i == (parttype+1)
        end # info.is_present[i] == Int32(1)
    end # i ∈ 1:size(info.is_present)[1]

    close(f)

end


"""
    read_block_data(f::IOStream, data_type::DataType, n_dim::Integer, npart::Integer)

Reads the binary data in a block.
"""
function read_block_data(f::IOStream, data_type::DataType, n_dim::Integer, npart::Integer)
    
    if n_dim > 1
        return copy(transpose(
        read!(f, Array{data_type,2}(undef, (n_dim,npart) ) 
                )))
        #return read!(f, Array{data_type,2}(undef, (npart,n_dim)))
        # return read!(f, Array{data_type,2}(undef, (n_dim,npart)))
    else
        read!(f, Array{data_type,1}(undef, npart))
    end
end

"""
    check_info(filename::String, blockname::String)

Helper function to read INFO block or construct `InfoLine` for MASS block, if no INFO block is present.
Returns a single `InfoLine` struct.
"""
function check_info(filename::String, blockname::String)
    info = read_info(filename)
    if info == 1
        if blockname == "MASS"
            return InfoLine("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
        else
            error("No Info block in snapshot! Supply InfoLine type!")
        end
    else # info != 1
        for i ∈ 1:size(info)[1]
            if info[i].block_name == blockname
                return info[i]
            end # if block found
        end # loop over info
        if isa(info, Array)
            if (blockname == "MASS")
                return InfoLine("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
            else
                error("Block not present!")
            end
        end
    end # info == 1
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
        return block_positions_dict[blockname]
    else
        if blockname != "MASS"
            # if the block is not present we need error handling!
            error("Requested block $blockname not present!")
        else
            return -1
        end

    end
end




"""
    read_block_with_offset(filename::String, data_old, pos0::Integer, info::InfoLine,
                                offset::Integer, offset_key, n_to_read::Integer, part_per_key )

Read part of a block and return an array of the entries.
"""
function read_block_with_offset(filename::String, data_old, pos0::Integer, info::InfoLine,
                                offset::Integer, offset_key, n_to_read::Integer, part_per_key )

    # open the file
    f = open(filename)

    # number of bits in data_type
    len = sizeof(info.data_type) * info.n_dim

    # jump to position of particle type in relevant block
    seek(f, pos0+offset*len)

    # store position in file
    p = position(f)

    # allocate array to store data
    data =  Array{info.data_type,2}(undef, (n_to_read, info.n_dim))

    n_read = 1
    n_this_key = 0

    for i = 1:size(offset_key)[1]

        # jump to start of key
        seek(f, p + len*offset_key[i])
        n_this_key += part_per_key[i]

        data[n_read:n_this_key, :] = read_block_data(f, info.data_type, info.n_dim, part_per_key[i])

        n_read += part_per_key[i]

    end # for

    # close the file
    close(f)

    return [data_old; data]
end

"""
    read_block_with_offset!(data, n_read::Integer, filename::String, pos0::Integer, info::InfoLine,
                                offset::Integer, offset_key::Array{<:Integer}, 
                                part_per_key::Array{<:Integer} )

Read part of a block to pre-allocated array.
"""
function read_block_with_offset!(data, n_read::Integer, filename::String, pos0::Integer, info::InfoLine,
                                offset::Integer, offset_key::Array{<:Integer}, 
                                part_per_key::Array{<:Integer} )

    # open the file
    f = open(filename)

    # number of bits in data_type
    len = sizeof(info.data_type) * info.n_dim

    # jump to position of particle type in relevant block
    seek(f, pos0+offset*len)

    # store position in file
    p = position(f)

    n_this_key = n_read 

    for i = 1:size(offset_key)[1]

        # jump to start of key
        seek(f, p + len*offset_key[i])
        n_this_key += part_per_key[i]

        data[n_read+1:n_this_key, :] = read_block_data(f, info.data_type, info.n_dim, part_per_key[i])

        n_read += part_per_key[i]

    end # for

    # close the file
    close(f)

    data
end