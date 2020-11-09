using Dates
using Base.Threads

"""
    struct KeyHeader

Helper struct to store header information of .key files.

# Fields

| Name                   | Meaning                                        | 
| :---                   | :------                                        |
| `nkeys_file`           | Number of keys in this .key file               |
| `domain_corners`       | Corners of the domains defined by these keys   |
| `domain_fac`           | Factor needed for reconstructung int positions |
| `bits`                 | Size of PH keys it bits                        |
| `nkeys_total`          | Total number of keys in all files              |
| `nkeys_total_highword` | Total number of keys in all files              |
"""
struct KeyHeader
    nkeys_file::Vector{Int32}
    domain_corners::Vector{Float64}
    domain_fac::Float64
    bits::Int32
    nkeys_total::Vector{UInt32}
    nkeys_total_highword::Vector{UInt32}
end

"""
    read_keyheader(filename::String)

Reads the header of a .key file.
"""
function read_keyheader(filename::String)

    f = open(filename)
    blocksize = read(f, Int32)

    if blocksize == 8
        swap = 0
        snap_format = 2
    elseif blocksize == 256
        swap = 0
        snap_format = 1
    else
        blocksize = bswap(blocksize)
        if blocksize == 8
            swap = 1
            snap_format = 2
        elseif blocksize == 256
            swap = 1
            snap_format = 1
        else
            error("incorrect file format encountered when reading header of $filename")
        end
    end

    if snap_format == 2
        seek(f, 16)
        skip_line = read(f, Int32)
    end

    nkeys = read!(f, Array{Int32,1}(undef,6))
    domain_corners = read!(f, Array{Float64,1}(undef,3))
    domain_fac = read(f, Float64)

    read!(f, Array{Float64,1}(undef,2))
    read(f, Float64)
    read(f, Float64)
    read(f, Int32)

    bits = read(f, Int32)
    nkeys_total = read!(f, Array{UInt32,1}(undef,6))

    read(f, Int32)
    read(f, Int32)
    read(f, Float64)
    read(f, Float64)
    read(f, Float64)
    read(f, Float64)
    read(f, Int32)
    read(f, Int32)

    nkeys_total_highword = read!(f, Array{UInt32,1}(undef,6))

    close(f)

    return KeyHeader( nkeys, domain_corners, domain_fac, bits, nkeys_total, nkeys_total_highword )
end


"""
    read_key_index(file_key_index::String)

Reads the .key.index file.
"""
function read_key_index(file_key_index::String)

    finfo = stat(file_key_index)

    f = open(file_key_index)
    n = read(f, Int32)

    if finfo.size == 4+n*8*2+n*4
        int_type = UInt64
    else
        int_type = UInt32
    end

    low_list = read!(f, Array{int_type,1}(undef,n))
    high_list = read!(f, Array{int_type,1}(undef,n))
    file_list = read!(f, Array{UInt32,1}(undef,n))

    close(f)

    return low_list, high_list, file_list
end