"""
    These functions read particles in a given range over multiple files, based
    on the relevant peano-hilbert keys.
    This code is based on read_particles_in_box.pro by Dr. habil. Klaus Dolag.
"""

using Dates
using Base.Threads

"""
    Functionality to read the header of the key files
"""

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
    Helper function for Peano Hilbert stuff
"""

"""
    get_int_pos(pos::Real, domain_corner::Real, domain_fac::Real )

Computes the integer position along the PH line.
"""
@inline function get_int_pos(pos::Real, domain_corner::Real, domain_fac::Real )
    return floor(Int64, ( pos - domain_corner ) * domain_fac) 
end

"""
    Lookup tables for peano hilbert keys
"""
const global rottable3 = [36  28  25  27  10  10  25  27;
                          29  11  24  24  37  11  26  26;
                           8   8  25  27  30  38  25  27;
                           9  39  24  24   9  31  26  26;
                          40  24  44  32  40   6  44   6;
                          25   7  33   7  41  41  45  45;
                           4  42   4  46  26  42  34  46;
                          43  43  47  47   5  27   5  35;
                          33  35  36  28  33  35   2   2;
                          32  32  29   3  34  34  37   3;
                          33  35   0   0  33  35  30  38;
                          32  32   1  39  34  34   1  31;
                          24  42  32  46  14  42  14  46;
                          43  43  47  47  25  15  33  15;
                          40  12  44  12  40  26  44  34;
                          13  27  13  35  41  41  45  45;
                          28  41  28  22  38  43  38  22;
                          42  40  23  23  29  39  29  39;
                          41  36  20  36  43  30  20  30;
                          37  31  37  31  42  40  21  21;
                          28  18  28  45  38  18  38  47;
                          19  19  46  44  29  39  29  39;
                          16  36  45  36  16  30  47  30;
                          37  31  37  31  17  17  46  44;
                          12   4   1   3  34  34   1   3;
                           5  35   0   0  13  35   2   2;
                          32  32   1   3   6  14   1   3;
                          33  15   0   0  33   7   2   2;
                          16   0  20   8  16  30  20  30;
                           1  31   9  31  17  17  21  21;
                          28  18  28  22   2  18  10  22;
                          19  19  23  23  29   3  29  11;
                           9  11  12   4   9  11  26  26;
                           8   8   5  27  10  10  13  27;
                           9  11  24  24   9  11   6  14;
                           8   8  25  15  10  10  25   7;
                           0  18   8  22  38  18  38  22;
                          19  19  23  23   1  39   9  39;
                          16  36  20  36  16   2  20  10;
                          37   3  37  11  17  17  21  21;
                           4  17   4  46  14  19  14  46;
                          18  16  47  47   5  15   5  15;
                          17  12  44  12  19  6   44   6;
                          13   7  13   7  18  16  45  45;
                           4  42   4  21  14  42  14  23;
                          43  43  22  20   5  15   5  15;
                          40  12  21  12  40   6  23   6;
                          13   7  13   7  41  41  22  20 ]

const global subpix3 = [0  7  1  6  3  4  2  5 ;
                        7  4  6  5  0  3  1  2 ;
                        4  3  5  2  7  0  6  1 ;
                        3  0  2  1  4  7  5  6 ;
                        1  0  6  7  2  3  5  4 ;
                        0  3  7  4  1  2  6  5 ;
                        3  2  4  5  0  1  7  6 ;
                        2  1  5  6  3  0  4  7 ;
                        6  1  7  0  5  2  4  3 ;
                        1  2  0  3  6  5  7  4 ;
                        2  5  3  4  1  6  0  7 ;
                        5  6  4  7  2  1  3  0 ;
                        7  6  0  1  4  5  3  2 ;
                        6  5  1  2  7  4  0  3 ;
                        5  4  2  3  6  7  1  0 ;
                        4  7  3  0  5  6  2  1 ;
                        6  7  5  4  1  0  2  3 ;
                        7  0  4  3  6  1  5  2 ;
                        0  1  3  2  7  6  4  5 ;
                        1  6  2  5  0  7  3  4 ;
                        2  3  1  0  5  4  6  7 ;
                        3  4  0  7  2  5  1  6 ;
                        4  5  7  6  3  2  0  1 ;
                        5  2  6  1  4  3  7  0 ;
                        7  0  6  1  4  3  5  2 ;
                        0  3  1  2  7  4  6  5 ;
                        3  4  2  5  0  7  1  6 ;
                        4  7  5  6  3  0  2  1 ;
                        6  7  1  0  5  4  2  3 ;
                        7  4  0  3  6  5  1  2 ;
                        4  5  3  2  7  6  0  1 ;
                        5  6  2  1  4  7  3  0 ;
                        1  6  0  7  2  5  3  4 ;
                        6  5  7  4  1  2  0  3 ;
                        5  2  4  3  6  1  7  0 ;
                        2  1  3  0  5  6  4  7 ;
                        0  1  7  6  3  2  4  5 ;
                        1  2  6  5  0  3  7  4 ;
                        2  3  5  4  1  0  6  7 ;
                        3  0  4  7  2  1  5  6 ;
                        1  0  2  3  6  7  5  4 ;
                        0  7  3  4  1  6  2  5 ;
                        7  6  4  5  0  1  3  2 ;
                        6  1  5  2  7  0  4  3 ;
                        5  4  6  7  2  3  1  0 ;
                        4  3  7  0  5  2  6  1 ;
                        3  2  0  1  4  5  7  6 ;
                        2  5  1  6  3  4  0  7 ]


"""
    peano_hilbert_key(bits::Integer, x::Integer, y::Integer, z::Integer)

Computes a Peano-Hilbert key for an integer triplet (x,y,z) with x,y,z in the range between 0 and 2^bits-1.
"""
function peano_hilbert_key(bits::Integer, x::Integer, y::Integer, z::Integer)

    rotation = 1
    key = 0
    mask = 1 << (bits-1)

    for imask = 1:bits

        pix = (
                (( x & mask) > 0 ? 4 : 0)  +
                (( y & mask) > 0 ? 2 : 0)  +
                (( z & mask) > 0 ? 1 : 0)
              ) + 1

        key = key << 3

        key = key | subpix3[rotation, pix]

        rotation = rottable3[rotation, pix] + 1

        mask = mask >> 1

    end

    return key
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

function get_index_bounds(ids::Vector{<:Integer}, low_bounds::Vector{<:Integer}, high_bounds::Vector{<:Integer})

    nids = length(ids)
    nbounds = length(low_bounds)

    ind_all = zeros(Int64, nids)

    icountall = 1

    icountids    = 1
    icountbounds = 1

    lend = false

    while !lend

        if low_bounds[icountbounds] <= ids[icountids] <= high_bounds[icountbounds]

            ind_all[icountall] = icountbounds
            icountall += 1

            lend2 = false
            while !lend2
                if ids[icountids] > high_bounds[icountbounds]
                    lend2 = true
                else # ids[icountids] > high_bounds[icountbounds]
                    icountids += 1

                    if icountids >= nids
                        lend2 = true
                    end # icountids >= nids
                end # ids[icountids] > high_bounds[icountbounds]
            end # while !lend2

            icountbounds += 1

        else # low_bounds[icountbounds] <= ids[icountids] <= high_bounds[icountbounds]

            if ids[icountids] < low_bounds[icountbounds]

                icountids += 1

                if icountids <= nids

                    lend2 = false
                    while !lend2

                        if ids[icountids] >= low_bounds[icountbounds]
                            lend2 = true
                        else # ids[icountids] >= low_bounds[icountbounds]
                            icountids += 1

                            if icountids >= nids
                                lend2 = true
                            end
                        end # if ids[icountids] >= low_bounds[icountbounds]

                    end # while !lend2

                else # if icountids < nids

                    if ids[icountids] > high_bounds[icountbounds]

                        icountbounds += 1

                        if icountbounds < nbounds

                            lend2 = false

                            while !lend2

                                if ids[icountids] <= high_bounds[icountbounds]
                                    lend2 = true
                                else # if ids[icountids] <= high_bounds[icountbounds]
                                    icountbounds += 1
                                    if icountbounds >= nbounds
                                        lend2 = true
                                    end
                                end # if ids[icountids] <= high_bounds[icountbounds]
                            end # while !lend2
                        end # if icountbounds < nbounds

                    end # if ids[icountids] > high_bounds[icountbounds]

                end # if icountids < nids

            else # if ids[icountids] < low_bounds[icountbounds]
                icountbounds += 1

                if icountbounds < nbounds
                    lend2 = false

                    while !lend2
                        if ids[icountids] <= high_bounds[icountbounds]
                            lend2 = true
                        else # ids[icountids] <= high_bounds[icountbounds]
                            icountbounds += 1
                            if icountbounds >= nbounds
                                lend2 = true
                            end # icountbounds >= nbounds
                        end # ids[icountids] <= high_bounds[icountbounds]
                    end # while !lend2

                end # icountbounds < nbounds

            end # if ids[icountids] < low_bounds[icountbounds]

        end # if low_bounds[icountbounds] <= ids[icountids] <= high_bounds[icountbounds]

        if icountids >= nids
            lend = true
        end
        if icountbounds >= nbounds
            lend = true
        end

        #
        # println("icountids     = ", icountids)
        # println("icountbounds  = ", icountbounds)

    end # while !lend

    if icountall > 1
        ind_out = ind_all[1:icountall]
        if ind_out[end] == 0
            return ind_out[1:end-1]
        else
            return ind_out
        end
    else
        return -1
    end
end

"""
    find_files_for_keys(filebase::String, nfiles::Integer, keylist::Vector{<:Integer})

Selects the files in which the particles associated with the given Peano-Hilbert keys are stored.
"""
function find_files_for_keys(filebase::String, nfiles::Integer, keylist::Vector{<:Integer})

    file_key_index = filebase * ".key.index"

    # if index file does not exist all key files need to be read
    if !isfile(file_key_index)
        return collect(0:nfiles-1)
    end

    # get the data from the index file
    low_list, high_list, file_list = read_key_index(file_key_index)

    # sort the keys
    key_sort = sortperm(keylist)

    index_bounds = get_index_bounds(keylist[key_sort], low_list, high_list)

    file_sort = sort(file_list[index_bounds])

    return Int64.(file_sort)
end

"""
    find_files_for_keys_AR(filebase::String, nfiles::Integer, keylist::Vector{<:Integer})

Selects the files in which the particles associated with the given Peano-Hilbert keys are stored. Version of Antonio. (Slower than Klaus' version!)
"""
function find_files_for_keys_AR(filebase::String, nfiles::Integer, keylist::Vector{<:Integer})

    file_key_index = filebase * ".key.index"

    # if index file does not exist all key files need to be read
    if !isfile(file_key_index)
        return collect(0:nfiles-1)
    end

    # get the data from the index file
    low_list, high_list, file_list = read_key_index(file_key_index)

    mask = falses(length(low_list))

    for key in keylist
        @. mask = mask | ( (key >= low_list ) & ( key <= high_list ))
    end

    return Int64.(unique!(sort!(file_list[mask])))
end


"""
    get_index_list(idarr1::Array{<:Integer}, idarr2::Array{<:Integer})

Get positions in `idarr2` where `idarr2` matches `idarr1`.
"""
@inline function get_index_list(idarr1::Array{<:Integer}, idarr2::Array{<:Integer})

    narr1 = length(idarr1)
    narr2 = length(idarr2)

    ind_all   = zeros(Int64, narr1)
    not_arr2t = zeros(Int64, narr1)
    not_arr1t = zeros(Int64, narr2)

    icountall     = 1
    icountarr1    = 1
    icountarr2    = 1
    icountnotarr1 = 1
    icountnotarr2 = 1

    lend = false

    iiarr2 = sortperm(idarr2[:,1])

    while !lend

        if idarr2[iiarr2[icountarr2]] == idarr1[icountarr1]

            ind_all[icountall] = iiarr2[icountarr2]
            icountall  += 1
            icountarr1 += 1
            icountarr2 += 1
        else  # idarr2[iiarr2[icountnotarr2]] == idarr1[icountnotar1]]
            if idarr2[iiarr2[icountarr2]] < idarr1[icountarr1]
                not_arr1t[icountnotarr1] = iiarr2[icountarr2]
                icountarr2    += 1
                icountnotarr1 += 1
            else # idarr2[iiarr2[icountnotarr2]] < idarr1[icountnotar1]]
                not_arr2t[icountnotarr2] = icountarr1
                icountarr1    += 1
                icountnotarr2 += 1
            end # idarr2[iiarr2[icountnotarr2]] < idarr1[icountnotar1]]
        end # idarr2[iiarr2[icountnotarr2]] == idarr1[icountnotar1]]
        if (icountarr2 >= narr2 ) || (icountarr1 >= narr1)
            lend = true
        end
    end # while !lend

    rest = narr1 - icountarr1
    if rest > 0
        not_arr2t[icountnotarr2:icountnotarr2+rest] = icountarr1:icountarr1+rest
        icountnotarr2 += rest
    end # rest > 0

    if icountall > 1
        ind_out = ind_all[1:icountall-1]
    # else
    #     error("Not enough data found!")
    end

    return ind_out
end


"""
    get_index_list_dict(keylist::Array{<:Integer}, keys_in_file::Array{<:Integer})

Get positions in `keys_in_file` where `keys_in_file` matches `keylist`. Uses a `Dict` for lookup -> slower than the normal version.
"""
@inline function get_index_list_dict(keylist::Array{<:Integer}, keys_in_file::Array{<:Integer})

    dict = Dict((n, i) for (i, n) in enumerate(keys_in_file))
    result = Vector{Int}(undef, length(keylist))
    len = 0

    for k in keylist
        i = get(dict, k, nothing)
        if i !== nothing
            len += 1
            @inbounds result[len] = i
        end
    end
    return resize!(result, len)
end

"""
    get_keylist(h_key::KeyHeader, x0::Array{<:Real}, x1::Array{<:Real})

Get all Peano-Hilbert keys for domain defined by the corner points `x0` and `x1`.
"""
function get_keylist(h_key::KeyHeader, x0::Array{<:Real}, x1::Array{<:Real})

    ix0 = zeros(Int, 3)
    ix1 = zeros(Int, 3)
    dix = zeros(Int, 3)

    nkeys = 1

    @inbounds for i = 1:3
        ix0[i] = get_int_pos( x0[i], h_key.domain_corners[i], h_key.domain_fac )
        ix1[i] = get_int_pos( x1[i], h_key.domain_corners[i], h_key.domain_fac )
        dix[i] = ix1[i] - ix0[i] + 1
        nkeys *= dix[i]
    end

    keylist = zeros(Int, nkeys)

    i = 1
    @inbounds for ix = ix0[1]:ix1[1], iy = ix0[2]:ix1[2], iz = ix0[3]:ix1[3]
        keylist[i] = peano_hilbert_key(h_key.bits, ix, iy, iz)
        i += 1
    end

    return sort(keylist)
end


"""
   function join_blocks(offset_key, part_per_key)
    
Joins neigboring blocks to simplify read-in.
"""
@inline function join_blocks(offset_key, part_per_key)

    use_block = trues(length(offset_key))

    icount = 1

    for i = 2:length(offset_key)-1

        if offset_key[i] == offset_key[icount] + part_per_key[icount]
            part_per_key[icount] += part_per_key[i]
            use_block[i] = false
        else
            icount += 1
        end # if

    end # for

    return use_block, part_per_key
end


"""
    find_read_positions(files::Array{<:Integer}, filebase::String, 
                             blocks::Array{String}, parttype::Integer)

Helper function to get positions and length of particle blocks in files.
"""
function find_read_positions(files::Array{<:Integer}, filebase::String, 
                             blocks::Array{String}, parttype::Integer,
                             keylist::Array{<:Integer}, key_info::Array{Info_Line},
                             verbose::Bool)

    # store number of file
    N_files = length(files)

    # allocate arrys to store reading information
    file_offset_key      = Array{Array{<:Integer}}(undef, N_files)
    file_part_per_key    = Array{Array{<:Integer}}(undef, N_files)
    file_block_positions = Array{Dict{String,Integer}}(undef, N_files)

    @inbounds for i = 1:N_files

        filename = select_file(filebase, files[i])

        # read header of the file
        h = head_to_obj(filename)

        if h.npart[parttype+1] == 0
            error("No particles of type $parttype in file!")
        end

        filename_keyfile = filename * ".key"

        # read key file data
        h_key = read_keyheader(filename_keyfile)
        keys_in_file = read_block_by_name(filename_keyfile, "KEY",
                                          info = key_info[getfield.(key_info, :block_name) .== "KEY"][1],
                                          parttype = parttype)

        if verbose
            @info "Calculating index list..."
            t1 = Dates.now()
        end

        index_list = get_index_list_dict(keylist, keys_in_file)

        if verbose
            t2 = Dates.now()
            @info "Index list done. Took: $(t2 - t1)"
            @info "Reading $(length(index_list)) key segments..."
        end

        # number of particles associated with PH key
        part_per_key = read_block_by_name(filename_keyfile, "NKEY",
                                          info = key_info[getfield.(key_info, :block_name) .== "NKEY"][1],
                                          parttype = parttype)


        # offsets in the blocks to get to the relevant particles
        offset_key = read_block_by_name(filename_keyfile, "OKEY",
                                          info = key_info[getfield.(key_info, :block_name) .== "OKEY"][1],
                                          parttype = parttype)

        # sort the offset arrays for simplification step
        sorted_offset = sortperm(offset_key[index_list])

        # save sorted arrays
        offset_key   = offset_key[index_list[sorted_offset]]
        part_per_key = part_per_key[index_list[sorted_offset]]

        # check if blocks can be joined
        use_block, part_per_key = join_blocks(offset_key, part_per_key)

        if verbose
            @info "Reduced independent blocks from $(length(offset_key)) to $(length(use_block[use_block]))"
        end

        # store the arrays for later reading
        file_offset_key[i]   = offset_key[use_block]
        file_part_per_key[i] = part_per_key[use_block]

        file_block_positions[i] = get_block_positions(filename)

        # double-check if blocks are present
        blocks_in_file = String.(keys(file_block_positions[i]))
        for blockname in blocks
            if !block_present(filename, blockname, blocks_in_file)
                error("Block $blockname not present in file $filename !")
            end
        end

    end # for i = 1:length(files)

    return file_offset_key, file_part_per_key, file_block_positions
end


"""
    Read particles in box main
"""



"""
    read_particles_in_box(filename::String, blocks::Vector{String},
                          corner_lowerleft, corner_upperright;
                          parttype::Integer=0, verbose::Bool=true)

Reads all particles within a box defined by a lower left and upper right corner
for a given particle type. Returns a dictionary with all requested blocks.
"""
function read_particles_in_box(filename::String, blocks::Vector{String},
                               corner_lowerleft, corner_upperright;
                               parttype::Integer=0, verbose::Bool=true)


    if verbose
        println()
        @info "Running on $(nthreads()) threads"
    end

    filebase = filename
    file_key_index = filebase * ".key.index"

    # if the snapshot does not exist it may be split into multiple files
    if !isfile(filebase)

        if verbose
            @info "File: $filebase not found, looking for sub-files."
        end

        # try reading the first of the distributed snapshots
        filename = select_file(filebase, 0)

        h = head_to_obj(filename)

        if verbose
            @info "$(h.num_files) sub-files found."
        end

        nfiles = h.num_files
    else
        nfiles = 1
        h = head_to_obj(filename)
    end

    blocks, no_mass_block = check_blocks(filename, blocks)

    # read info blocks once here
    snap_info = read_info(filename)
    key_info  = read_info(filename * ".key")

    if verbose
        @info "All requested blocks present!"
        @info "Checking for .key files..."
    end

    # check if key files are present
    file_key = filename * ".key"
    if !block_present(file_key, "KEY")
        error("No .key file present!")
    end

    if verbose
        @info ".key files found!"
        @info "Calculating peano-hilbert keys..."
        t1 = Dates.now()
    end

    # first read the header
    h_key = read_keyheader(file_key)

    # get a list of the required peano-hilbert keys
    keylist = get_keylist(h_key, corner_lowerleft, corner_upperright)

    if verbose
        t2 = Dates.now()
        @info "$(length(keylist)) Peano-Hilbert keys found. Took: $(t2 - t1)"
        @info "Looking for relevant files..."
        t1 = Dates.now()
    end

    # find relevant files
    files = find_files_for_keys(filebase, nfiles, keylist)
    
    N_files = length(files)
    if verbose
        t2 = Dates.now()
        @info "$N_files files found. Took: $(t2 - t1)"

        @info "Searching read positions..."
        println()
        t1 = Dates.now()
    end

    # find all the positions where to read data
    file_offset_key, file_part_per_key, file_block_positions = find_read_positions( files, filebase, blocks, 
                                                                                    parttype, keylist, key_info, 
                                                                                    verbose)

    N_to_read = 0

    @inbounds for i = 1:N_files
        N_to_read += sum(file_part_per_key[i])
    end

    if verbose
        t2 = Dates.now()
        println()
        @info "Positions read. Took: $(t2 - t1)"
        println()
        @info "Reading $N_to_read particles..."
    end

    # prepare dictionary for particle storage
    d = allocate_data_dict(blocks, N_to_read, snap_info, no_mass_block)

    if verbose
        @info "Reading Blocks..."
        t1 = Dates.now()
    end

    n_read = 1

    for i = 1:N_files
        if N_files > 1
            filename = filebase * ".$(files[i])"
        else
            filename = filebase
        end

        # read blocks in parallel
        @threads for j = 1:length(blocks)

            block_info = snap_info[getfield.(snap_info, :block_name) .== blocks[j]][1]

            # add offset of particle types that should not be read
            offset = 0
            for i=1:length(h.npart)
                if (block_info.is_present[i] > 0) & (h.npart[i] > 0) & ( i < parttype + 1)
                    offset += h.npart[i]
                end
            end

            # reads data into the dictionary and counts up n_read
            read_block_with_offset!(d[blocks[j]], n_read, filename, 
                                    file_block_positions[i][blocks[j]],
                                    block_info, offset, file_offset_key[i],
                                    file_part_per_key[i])

        end # loop over blocks

        
        n_read += sum(file_part_per_key[i])

        @info "Read $(n_read-1) / $N_to_read particles"

    end # for i = 1:length(files)

    # finally construct masses of no mass block present
    if no_mass_block
        d["MASS"] = h.massarr[parttype+1] .* ones(Float32, N_to_read, 1)
    end

    if verbose
        t2 = Dates.now()
        @info "Blocks read. Took: $(t2 - t1)"
    end

    return d
end


"""
    read_particles_in_box(filename::String, blocks::String,
                          corner_lowerleft, corner_upperright;
                          parttype::Integer=0, verbose::Bool=true)

Like `read_particles_in_box` but for a single block. Returns the block as an array.

See also: [`read_particles_in_box`](@ref)
"""
function read_particles_in_box(filename::String, blocks::String,
                               corner_lowerleft, corner_upperright;
                               parttype::Integer=0, verbose::Bool=true)

    d = read_particles_in_box(filename, [blocks], corner_lowerleft, corner_upperright, parttype=parttype, verbose=verbose)

    return d[blocks]
end


"""
    read_particles_in_box(filename::String, blocks::Vector{String},
                          center_pos, radius;
                          parttype::Integer=0, verbose::Bool=true)

Reads all particles within a box encapsulating a volume defined by center position
and radius for a given particle type. Returns a dictionary with all requested blocks.

See also: [`read_particles_in_box`](@ref)
"""
function read_particles_in_volume(filename::String, blocks::Vector{String},
                                  center_pos::Array{<:Real}, radius::Real;
                                  parttype::Integer=0, verbose::Bool=true)

    # calculate lower left and upper right corner
    x0 = center_pos .- radius
    x1 = center_pos .+ radius

    return read_particles_in_box(filename, blocks, x0, x1, parttype=parttype, verbose=verbose)
end

"""
    read_particles_in_box(filename::String, blocks::String,
                          center_pos, radius;
                          parttype::Integer=0, verbose::Bool=true)

Like `read_particles_in_volume` but for a single block. Returns the block as an array.

See also: [`read_particles_in_volume`](@ref)
"""
function read_particles_in_volume(filename::String, blocks::String,
                                  center_pos, radius;
                                  parttype::Integer=0, verbose::Bool=true)

    d = read_particles_in_volume(filename, [blocks], center_pos, parttype,
                                 parttype=parttype, verbose=verbose)

    return d[blocks]
end

