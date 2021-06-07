"""
            Functions to write snapshots.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

function write_format2_block_header(f::IOStream, blockname::String)

    # fill blockname to 4 characters
    if length(blockname) < 4
        padding = 4 - length(blockname)
        for i=1:padding
            blockname = blockname * " "
        end
    end

    # transform the name into an array of charcters.
    # Since Julia uses 2 byte chars and C uses 1 byte I have to take
    # a different approach and use Int8s.
    write_name =  Int8.(collect(blockname))

    @info "Writing block: $blockname"

    # write blocksite (8 bytes)
    write(f, UInt32(8))

    # write block name
    write(f, write_name)

    return f
end


"""
    write_block(f::IOStream, data,
                blockname::String="";
                snap_format::Integer=2)

Write `data` to a block to an opened file `f`.
"""
function write_block(f::IOStream, data,
                     blockname::String="";
                     snap_format::Integer=2)

    # get total number of particles to write
    N = size(data, ndims(data))

    # write blocksize
    dtype = typeof(data[1,1])
    dims = ndims(data)
    if dims == 2
        dims = size(data,1)
    end

    blocksize = UInt32(N * sizeof(dtype) * dims)

    if snap_format == 2

        # break if blockname not specified
        if blockname == ""
            error("Please specify blockname!")
        end

        # write header
        f = write_format2_block_header(f, blockname)

        # write end if name header, namely:
        # size to skip from directly after the block name and size of
        # the format 2 block header
        write(f, UInt32(blocksize + UInt32(8)))
        write(f, UInt32(8))
    end

    # write blocksize
    write(f, blocksize)

    # write the block. Since Julia stores the arrays differently in memory
    # they have to be transposed before the can be written.
    write(f, data)

    @info "Writing block done."

    write(f, blocksize)
end


"""
    write_header(f::IOStream, h::SnapshotHeader; snap_format::Integer=2)

Writes the header block to an opened file `f`.
"""
function write_header(f::IOStream, h::SnapshotHeader; snap_format::Integer=2)

    if snap_format == 2
        f = write_format2_block_header(f, "HEAD")
        write(f, UInt32(264))
        write(f, UInt32(8))
    end

    # write blocksize
    write(f, UInt32(256))

    # write header to file
    for fields in fieldnames(SnapshotHeader)
        write(f, getfield(h, fields))
    end

    # write blocksize
    write(f, UInt32(256))
end

"""
    convert_string_to_int8(name::String)

Converts a string to an array of Int8 to write to file.
"""
convert_string_to_int8(name::String) = [Int8(Char(name[i])) for i = 1:length(name)]

"""
    function get_datatype_name(dt::DataType, dim::Integer)

Convert the datatype into a string and then to a vector of 8bit integers to match name in Gadget.
"""
function get_datatype_name(dt::DataType, dim::Int32)

    # convert datatype into string
    if dt == Float32
        name = "FLOAT"
    elseif dt == Float64
        name = "DOUBLE"
    elseif dt == UInt32
        name = "LONG"
    elseif dt == UInt64
        name = "LLONG"
    end

    # add N if 3D
    if dim == 3
        name *= "N"
    end

    # pad to 8 characters
    for _ = length(name):7
        name *= " "
    end

    # return the datatype string as 8bit integers
    return convert_string_to_int8(name)
end

"""
    get_int8_blockname(info::InfoLine)

Pad the stored block name with spaces and convert to Array of Int8.
"""
function get_int8_blockname(info_line::InfoLine)

    block_name = info_line.block_name
        
    # pad block name with spaces
    while length(block_name) < 4
        block_name *= " "
    end

    return convert_string_to_int8(block_name)
end

"""
    write_info_block(f::IOStream, info::Vector{InfoLine}; snap_format::Integer=2)

Writes the info block to an opened file `f`.
"""
function write_info_block(f::IOStream, info::Vector{InfoLine}; snap_format::Integer=2)

    # every info line takes up 40 bits
    blocksize = UInt32(40 * length(info) )

    if snap_format == 2
        f = write_format2_block_header(f, "INFO")
        write(f, UInt32(blocksize+8))
        write(f, UInt32(8))
    end

    write(f, UInt32(blocksize))

    for i = 1:length(info)
        
        # get the block name in the correct format
        block_name = get_int8_blockname( info[i] )

        # get the datatype name in the correct format
        datatype_name = get_datatype_name( info[i].data_type, info[i].n_dim )

        # write the data
        write(f, block_name)
        write(f, datatype_name)
        write(f, info[i].n_dim)
        write(f, info[i].is_present)

    end

    write(f, UInt32(blocksize))
end
