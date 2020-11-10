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
                            info::InfoLine=InfoLine(),
                            parttype::Integer=-1)


    # read header - super fast and needed for flexibility
    h = head_to_obj(filename)

    if parttype == -1
        parttypes = ["PartType0", "PartType1", "PartType2",
                     "PartType3", "PartType4", "PartType5"]
        d = Dict()
    end # if parttype

    blockname = strip(blockname)

    if info.block_name == ""
        info = read_info(filename)
        if info == 1
            if blockname == "MASS"
                info = InfoLine("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
            else
                error("No Info block in snapshot! Supply InfoLine type!")
            end
        else # info != 1
            for i ∈ 1:size(info)[1]
                if info[i].block_name == blockname
                    info = info[i]
                    break
                end # if block found
            end # loop over info
            if isa(info, Array)
                if (blockname == "MASS")
                    info = InfoLine("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
                else
                    error("Block not present!")
                end
            end
        end # info == 1
    end

    f = open(filename)

    first_block = read(f, UInt32)

    if first_block != UInt32(8)
        error("Only snapshots of format 2 can be read by name!")
        return
    end

    while eof(f) != true
        blname = read_bockname(f)

        if blname == blockname
            break
        end
        p = position(f)
        seek(f,p+8)
        blocksize = read(f, UInt32)

        # check if the blocksize is correct
        blocksize = check_blocksize(f, p, blocksize)

        seek(f,p+blocksize+20)
    end

    if eof(f) != true
        p = position(f)
        seek(f,p+12)

        for i ∈ 1:size(info.is_present)[1]
            p = position(f)

            if info.is_present[i] == Int32(1)
                if parttype == -1
                    d[parttypes[i]] = Dict()
                    d[parttypes[i]][blockname] = read_block_data(f, info.data_type, info.n_dim, h.npart[i])
                else
                    if i == (parttype+1)

                        block = read_block_data(f, info.data_type, info.n_dim, h.npart[i])
                        close(f)
                        return block
                    else
                        seek(f, p + ( sizeof(info.data_type)*info.n_dim*h.npart[i] ))
                    end # if i == (parttype+1)
                end # parttype == -1
            end # info.is_present[i] == Int32(1)
        end # i ∈ 1:size(info.is_present)[1]

        # fill in mass from header
        if blockname == "MASS"
            for i ∈ 1:size(info.is_present)[1]
                if info.is_present[i] == Int32(0)
                    if parttype == -1
                        d[parttypes[i]] = Dict()
                        d[parttypes[i]][blockname] = Array{info.data_type,2}(undef,(h.npart[i], info.n_dim))
                        d[parttypes[i]][blockname] .= h.massarr[i]
                    end # parttype == -1
                end # info block
            end # loop over i
        end
        close(f)

        # this needs better error handling!
        return d

    else # blockname not found
        if blockname != "MASS"
            error("Block not present!")
        else
            if parttype == -1
                for i = 1:size(parttypes)[1]
                    d[parttypes[i]] = Dict()
                    d[parttypes[i]][blockname] = Array{info.data_type,2}(undef,(h.npart[i], info.n_dim))
                    d[parttypes[i]][blockname] .= h.massarr[i]
                end
                
                return d
            else
                block = Array{info.data_type,2}(undef,(h.npart[parttype+1], info.n_dim))
                block .= h.massarr[parttype+1]
                return block
            end # parttype == -1
        end
    end # eof(f) != true

end

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

function snap_2_d_info(filename::String, d::Dict{Any,Any}, info::Array{InfoLine,1})

    f = open(filename)

    seek(f, 300)

    for i ∈ 1:size(d["Header"]["npart"])[1]
        if d["Header"]["npart"][i] != Int32(0)
            d[d["Header"]["PartTypes"][i]] = Dict()
            d[d["Header"]["PartTypes"][i]]["MASS"] = Float32.(d["Header"]["massarr"][i] .* ones(d["Header"]["npart"][i]))
        end
    end

    for i ∈ 1:size(info)[1]
        #println("Reading block: ", info[i].block_name)
        d = read_block_with_info(f, d, info[i])
        p = position(f)
        seek(f, p+24)
    end

    close(f)

    return d

end

function snap_2_d(filename::String, data::Dict{Any,Any})

    f = open(filename)

    seek(f, 296)

    N = sum(data["Header"]["npart"])
    blocksize = read(f, UInt32)
    bit_size = Int64(blocksize/(3*N))

    if Int(bit_size) == 4
        dtype = Float32
        println("Reading single precision snapshot")
    elseif Int(bit_size) == 8
        dtype = Float64
        println("Reading double precision snapshot")
    else
        println("read error! neither 32 nor 64 bits data!")
        return -1
    end

    seek(f, 288)

    # set up dictionaries for particles
    for i in 1:size(data["Header"]["PartTypes"])[1]
        if data["Header"]["npart"][i] != 0
            data[data["Header"]["PartTypes"][i]] = Dict()
        end
        if data["Header"]["massarr"][i] != 0
           data[data["Header"]["PartTypes"][i]]["MASS"] = dtype.(data["Header"]["massarr"][i] * ones(data["Header"]["npart"][i]))
        end
    end

    blockname = "POS"


    while (blockname != "INFO") #& (eof(f) != true)

        #println(blockname)
        #println(typeof(blockname))
        # skip identifiers
        p = position(f)

        seek(f, p+8)

        if eof(f) == true
            break
        end

        p, data = read_block(p, data, dtype, blockname, f, bit_size)


        seek(f,p+4)

        if eof(f) == true
            break
        end

        p = position(f)
        seek(f, p+4)
        # read blockname

        blockname = read_bockname(f)

        p = position(f)
        #seek(f, p+8)

    end


    close(f)

    return data


end

function read_block_with_info(f::IOStream, d::Dict{Any,Any}, info::InfoLine)

    parttypes = ["PartType0", "PartType1", "PartType2",
                 "PartType3", "PartType4", "PartType5"]

    for i ∈ 1:size(info.is_present)[1]

        if info.is_present[i] == Int32(1)
            d[parttypes[i]][info.block_name] = read_block_data(f, info.data_type, info.n_dim, d["Header"]["npart"][i])
        end

    end

    return d

end

function read_block(p::Integer, data::Dict{Any,Any}, dtype::DataType, blockname::String,
                    f::IOStream, bit_size::Integer)



    file_curr = @__FILE__
    path_curr = file_curr[1:end-17]

    include(path_curr * "part_specific_fields.jl")

    N = sum(data["Header"]["npart"])

    blocksize = read(f, UInt32)

    if blockname == "MASS"
        for i ∈ 1:size(data["Header"]["PartTypes"])[1]

            n = data["Header"]["npart"][i]

            if data["Header"]["npart"][i] != Int32(0)

                if data["Header"]["massarr"][i] == Int32(0)
                    data[data["Header"]["PartTypes"][i]][blockname] = read_block_data(f, dtype, 1, n)
                else
                    data[data["Header"]["PartTypes"][i]][blockname] = dtype.(data["Header"]["massarr"][i] .* ones(n))
                end

            end
        end
    elseif blockname != "ID"

        gas_block = any(x->x==blockname,gas_arr)

        if gas_block == true

            #println("Reading gas block")

            n = Int64(data["Header"]["npart"][1])

            dim = Int(blocksize/(n*bit_size))

            #println("nr. of gas-particles: ", n)
            #println("reading dimensions: ", blocksize/(n*bit_size))

            data[data["Header"]["PartTypes"][1]][blockname] = read_block_data(f, dtype, dim, n)

        # check if bh
        elseif (any(x->x==blockname,bh_arr)) == true

            n = Int64(data["Header"]["npart"][6])

            #println("nr. of bh-particles: ", n)

            dim = Int(blocksize/(n*bit_size))

            #println("reading dimensions: ", blocksize/(n*bit_size))

            if blockname != "BHPC"
                data[data["Header"]["PartTypes"][6]][blockname] = read_block_data(f, dtype, dim, n)
            else
                data[data["Header"]["PartTypes"][6]][blockname] = Int64.(read_block_data(f, UInt32, 1, n))
            end

        elseif blockname == gas_bh

            n = Int64(data["Header"]["npart"][1])
            data[data["Header"]["PartTypes"][1]][blockname] = Int64.(read_block_data(f, UInt32, 1, n))

            n = Int64(data["Header"]["npart"][6])
            data[data["Header"]["PartTypes"][6]][blockname] = Int64.(read_block_data(f, UInt32, 1, n))

        else
            # read Blocks
            #println("Blockname: ", blockname)

            for i in 1:size(data["Header"]["PartTypes"])[1]

                dim = (blocksize/(N*bit_size))

                dim = Int64(trunc(dim))

                if data["Header"]["npart"][i] != 0

                    n = Int64(data["Header"]["npart"][i])

                    data[data["Header"]["PartTypes"][i]][blockname] = read_block_data(f, dtype, dim, n)

                end

            end
        end

    else

        for i in 1:size(data["Header"]["PartTypes"])[1]

            if data["Header"]["npart"][i] != 0

                n = Int64(data["Header"]["npart"][i])

                data[data["Header"]["PartTypes"][i]][blockname] = read_block_data(f, UInt32, 1, n)

            end

        end

    end

    p = position(f)

    return p, data
end


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