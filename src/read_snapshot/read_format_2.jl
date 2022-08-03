"""
    read_block(filename::String, blockname::String;
                                info::InfoLine=InfoLine(),
                                parttype::Integer=-1)

Reads a block in a snapshot with given name. Names are case sensitive.

# Examples
```julia
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
                    h::Union{Nothing,SnapshotHeader}=nothing,
                    return_haloid::Bool=false,
                    thisfilenum::Integer=0,
                    offset=0, nread=-1)

    # the file is actually a filebase -> return concatenated arrays
    if !isfile(filename)
        return read_block_subsnaps(filename, blockname; parttype, info, h, return_haloid)
    end

    # we need to know which particle to read
    if parttype == -1
        error("Please specify particle type!")
    end

    if isnothing(h)
        # read header - super fast and needed for flexibility
        h = read_header(filename)
    end


    # check if particle type is present
    if iszero(get_total_particles(h, parttype))
        error("Particle Type $parttype not present in simulation!")
    end

    # check if info is present
    if isnothing(info)
        info = check_info(filename, blockname)
    end

    if nread == -1
        nread = h.npart[parttype+1]
    end

    # check if a block position is supplied
    if block_position == -1

        block_position = check_block_position(filename, blockname)

        if block_position == -1 || (info.is_present[parttype+1] == 0)
            # if no mass block is present we can read it from the header
            if blockname == "MASS"
                block = Array{info.data_type,1}(undef, nread)
                block .= h.massarr[parttype+1]
                return block
            else
                error("Block $blockname not present for particle type $parttype !")
            end
        end
    end

    # allocate and fill HaloID array
    if return_haloid
        haloids = Vector{HaloID}(undef, h.npart[parttype+1])

        for j = 1:h.npart[parttype+1]
            haloids[j] = HaloID(thisfilenum, j)
        end
    end


    f = open(filename)

    seek(f,block_position)

    for i ∈ 1:size(info.is_present,1)
        p = position(f)

        if info.is_present[i] == Int32(1)
            if i == (parttype+1)
                skip(f, sizeof(info.data_type)*info.n_dim*offset)
                block = read_block_data(f, info.data_type, info.n_dim, nread)
                close(f)
                if return_haloid
                    return block, haloids
                else
                    return block
                end
            else
                seek(f, p + ( sizeof(info.data_type)*info.n_dim*h.npart[i] ))
            end # if i == (parttype+1)
        end # info.is_present[i] == Int32(1)
    end # i ∈ 1:size(info.is_present,1)

    close(f)

end

function read_block!(a::AbstractArray, filename::String, blockname::String;
                    parttype::Integer=-1,
                    block_position::Integer=-1,
                    info::Union{Nothing,InfoLine}=nothing,
                    h::Union{Nothing,SnapshotHeader}=nothing,
                    offset=0)
    # check if a block position is supplied
    if block_position == -1

        block_position = check_block_position(filename, blockname)

        if block_position == -1 || (info.is_present[parttype+1] == 0)
            # if no mass block is present we can read it from the header
            if blockname == "MASS"
                block = Vector{info.data_type}(undef, length(a))
                block .= h.massarr[parttype+1]
                return block
            else
                error("Block $blockname not present for particle type $parttype !")
            end
        end
    end


    f = open(filename)

    seek(f,block_position)

    for i ∈ 1:size(info.is_present,1)
        p = position(f)

        if info.is_present[i] == 1
            if i == (parttype+1)
                skip(f, sizeof(info.data_type) * info.n_dim * offset)
                read_block_data!(a, f, info.data_type, info.n_dim, -1) # need to remove excess parameters
                close(f)
                return a
            else
                seek(f, p + ( sizeof(info.data_type)*info.n_dim*h.npart[i] ))
            end # if i == (parttype+1)
        end # info.is_present[i] == Int32(1)
    end # i ∈ 1:size(info.is_present,1)

    close(f)

end

function read_block_prefiltered(filename::String, blockname::String, matched;
                    parttype::Integer=-1,
                    block_position::Integer=-1,
                    info::Union{Nothing,InfoLine}=nothing,
                    h::Union{Nothing,SnapshotHeader}=nothing)

    # the file is actually a filebase -> return concatenated arrays
    if !isfile(filename)
        return read_block_subsnaps(filename, blockname; parttype, info, h)
    end

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

        if block_position == -1 || (info.is_present[parttype+1] == 0)
            # if no mass block is present we can read it from the header
            if blockname == "MASS"
                block = Array{info.data_type,1}(undef, length(matched))
                block .= h.massarr[parttype+1]
                return block
            else
                error("Block $blockname not present for particle type $parttype !")
            end
        end
    end


    f = open(filename)

    seek(f,block_position)

    for i ∈ 1:size(info.is_present,1)
        p = position(f)

        if info.is_present[i] == Int32(1)
            if i == (parttype+1)
                block = read_block_data_prefiltered(f, info.data_type, info.n_dim, matched)
                close(f)
                return block
            else
                seek(f, p + ( sizeof(info.data_type)*info.n_dim*h.npart[i] ))
            end # if i == (parttype+1)
        end # info.is_present[i] == Int32(1)
    end # i ∈ 1:size(info.is_present,1)

    close(f)

end

"""
    read_block_subsnaps(filename::String, blockname::String;
                                info::InfoLine=InfoLine(),
                                parttype::Integer=-1)

Reads a block over all sub-snapshots. Names are case sensitive.

# Examples
```julia
julia> pos_info = InfoLine("POS", Float32, 1, [1, 1, 1, 1, 1, 1])
[...]
julia> gas_pos = read_block(filename, "POS", info=pos_info, parttype=0)
[...]
```
"""
function read_block_subsnaps(filebase::String, blockname::String;
                             parttype::Integer=-1,
                             info::Union{Nothing,InfoLine}=nothing,
                             h::Union{Nothing,SnapshotHeader}=nothing,
                             return_haloid::Bool=false)

    # the file is actually a filebase -> return concatenated arrays
    if !isfile(filebase * ".0")
        error("Neither $filebase nor $(filebase * ".0") present!")
    end

    # we need to know which particle to read
    if parttype == -1
        error("Please specify particle type!")
    end

    if isnothing(h)
        # read header - super fast and needed for flexibility
        h = read_header(filebase * ".0")
    end

    # check if info is present
    if isnothing(info)
        info = check_info(filebase * ".0", blockname)
    end

    # check if particle type is present
    if h.nall[parttype+1] == UInt32(0)
        error("Particle Type $parttype not present in simulation!")
    end

    return read_subsnaps(filebase, blockname, parttype, info, h; return_haloid)

end

"""
    read_subsnaps(filebase::String, blockname::String, parttype::Integer,
                          info::InfoLine, h_global::SnapshotHeader)

Reads a block over distributed files and returns it in one large Array.
"""
function read_subsnaps(filebase::String, blockname::String, parttype::Integer,
                          info::InfoLine, h_global::SnapshotHeader; return_haloid::Bool=false)

    # allocate empty array for block
    if info.n_dim > 1
        block = Array{info.data_type,2}(undef, ( info.n_dim, h_global.nall[parttype+1] ))
    else
        block = Array{info.data_type,1}(undef, h_global.nall[parttype+1])
    end

    # allocate HaloID array
    if return_haloid
        haloids = Vector{HaloID}(undef, h_global.nall[parttype+1])
    end

    # store number of particles that have been read
    N_read = 0

    # loop over all sub-files
    @inbounds for i = 0:(h_global.num_files-1)
        
        # get filename of the current file
        filename = "$filebase.$i"
        
        # read header of subfile
        h = read_header(filename)

        info = check_info(filebase * ".$i", blockname)
        if info.is_present[parttype+1] == 0
            continue
        end

        # get number of particles for this sub-file
        N_to_read = h.npart[parttype+1]

        # read the block (info has to be read individually for each file)
        if info.n_dim > 1
            read_block!(@view(block[:, (N_read+1):(N_read+N_to_read)]), filename, blockname; parttype, info, h)
        else
            read_block!(@view(block[(N_read+1):(N_read+N_to_read)]), filename, blockname; parttype, info, h)
        end

        # fill haloids array
        if return_haloid
            for j = 1:N_to_read
                haloids[N_read + j] = HaloID(i, j)
            end
        end

        # update number of read particles
        N_read += N_to_read
    end

    if return_haloid
        return block, haloids
    else
        return block
    end
end



"""
    read_block_data(f::IOStream, data_type::DataType, n_dim::Integer, npart::Integer)

Reads the binary data in a block.
"""
function read_block_data(f::IOStream, data_type::DataType, n_dim::Integer, npart::Integer)
    
    if n_dim > 1
        return read!(f, Array{data_type,2}(undef, (n_dim, npart)))
    else
        read!(f, Array{data_type,1}(undef, npart))
    end
end

function read_block_data!(a::AbstractArray, f::IOStream, data_type::DataType, n_dim::Integer, npart::Integer)
    read!(f, a)
end

function read_block_data_prefiltered(f::IOStream, data_type::DataType, n_dim::Integer, matched)
    firstind = first(matched)
    lastind = last(matched)

    if n_dim > 1
        a = Matrix{data_type}(undef, n_dim, lastind - firstind + 1)
    else
        a = Vector{data_type}(undef, lastind - firstind + 1)
    end

    skip(f, sizeof(data_type) * n_dim * (firstind - 1))

    return read!(f, a)
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
        for i ∈ 1:size(info,1)
            if info[i].block_name == blockname
                return info[i]
            end # if block found
        end # loop over info
        if isa(info, Array)
            if (blockname == "MASS")
                return InfoLine("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
            else
                error("Block $blockname not present!")
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

    for i = 1:size(offset_key,1)

        # jump to start of key
        seek(f, p + len*offset_key[i])
        n_this_key += part_per_key[i]

        # note the Int(...) are necessary to enable the reading into the SubArray returned by @view
        if info.n_dim == 1
            read_block_data!(@view(data[Int(n_read+1):Int(n_this_key)]), f, info.data_type, info.n_dim, part_per_key[i])
        else
            read_block_data!(@view(data[:, Int(n_read+1):Int(n_this_key)]), f, info.data_type, info.n_dim, part_per_key[i])
        end

        n_read += part_per_key[i]

    end # for

    # close the file
    close(f)

    data
end
