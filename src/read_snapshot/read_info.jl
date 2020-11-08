"""
    struct Info_Line([  block_name="", data_type=Float32, n_dim=Int32(0),
                                is_present=zeros(Int32, 6) ])

Contains the data of a single entry in the `INFO` block of a Gadget snapshot.

# Fields
| Name                                 | Meaning                                                                                |
| :----------------------------------  | :------------------------------------------------------------------------------------- |
| `block_name::String`                 | name of the data block, e.g. "POS"                                                     |
| `data_type::DataType`                | datatype of the block, e.g. Float32 for single precision, Float64 for double           |
| `n_dim::Int32`                       | number of dimensions of the block, usually 1 or 3                                      |
| `is_present::Vector{Int32}`          | array of flags for which particle type this block is present,                          |
|                                      |  e.g. gas only:  [ 1, 0, 0, 0, 0, 0 ],                                                 |
|                                      |  or gas + BHs: [ 1, 0, 0, 0, 0, 1 ]                                                    |

"""
struct Info_Line
    block_name::String              # name of the data block, e.g. "POS"
    data_type::DataType             # datatype of the block, e.g. Float32 for single precision, Float64 for double
    n_dim::Int32                    # number of dimensions of the block, usually 1 or 3
    is_present::Vector{Int32}       # array of flags for which particle type this block is present,
                                    # e.g. gas only:  [ 1, 0, 0, 0, 0, 0 ]
                                    # e.g. gas + BHs: [ 1, 0, 0, 0, 0, 1 ]

    function Info_Line(block_name="", data_type=Float32, n_dim=Int32(0),
                        is_present=zeros(Int32, 6))

        new(block_name, data_type, n_dim, is_present)
    end
end



"""
    read_info(filename; verbose::Bool=false)

Reads the info block of a snapshot and returns the information in an array of `Info_Line` types.
If `verbose=true` the blocknames are also printed to console.
"""
function read_info(filename::String; verbose::Bool=false)

    f = open(filename)
    blocksize = read(f, Int32)

    while eof(f) != true

        # read block name
        blockname = read_bockname(f)

        p = position(f)

        seek(f,p+8)

        skipsize = read(f, UInt32)

        skipsize = check_blocksize(f, p, skipsize)

        if blockname == "INFO"
            seek(f,p+12)
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
            seek(f,p+skipsize+20)
        end # if blockname == "INFO"

    end # while eof(f) != true

    println("No info block present!")

    return 1
end

"""
    read_info_line(f::IOStream)

Helper function to read the binary data into a `Info_Line` struct.
"""
function read_info_line(f::IOStream)

    block_name = read_bockname(f)

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