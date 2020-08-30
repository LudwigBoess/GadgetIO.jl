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
    write(f, Int32(8))

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


    # turn array of arrays d into one large array
    #data = reduce(vcat, d)     # OLD! Keep for later, maybe

    # get total number of particles to write
    N = length(data[:,1])

    # write blocksize
    dtype = typeof(data[1,1])
    dims = length(data[1,:])
    blocksize = Int32(N * sizeof(dtype) * dims)

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
        write(f, Int32(blocksize + Int32(8)))
        write(f, Int32(8))
    end

    # write blocksize
    write(f, blocksize)

    # write the block. Since Julia stores the arrays differently in memory
    # they have to be transposed before the can be written.
    write(f, copy(transpose(data)))

    @info "Writing block done."

    write(f, blocksize)
end


"""
    write_header(f::IOStream, h::Header; snap_format::Integer=2)

Writes the header block to an opened file `f`.
"""
function write_header(f::IOStream, h::Header; snap_format::Integer=2)

    if snap_format == 2
        f = write_format2_block_header(f, "HEAD")
        write(f, Int32(264))
        write(f, Int32(8))
    end

    # write blocksize
    write(f, Int32(256))

    # write header to file
    for fields in fieldnames(Header)
        write(f, getfield(h, fields))
    end

    # write blocksize
    write(f, Int32(256))
end