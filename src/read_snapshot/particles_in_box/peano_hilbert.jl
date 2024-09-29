"""
    get_int_pos(pos::Real, domain_corner::Real, domain_fac::Real, bits::Integer)

Computes the integer position along the PH line.
"""
@inline function get_int_pos(pos::Real, domain_corner::Real, domain_fac::Real, bits::Integer )
    val = (pos - domain_corner) * domain_fac 

    # round down, but negative values have to be reduced by one and values beyond the box border increased by one (unclear why exactly, probably due to how the bits work out, but this was carefully tested to be sure that this version works)
    val_int = floor(Int64, val) + ifelse(pos ≥ 0, 0, -1) + ifelse(val ≥ 2^bits - 0.5, 1, 0) 

    return val_int
end

# Lookup tables for peano hilbert keys
const rottable3 = [36  28  25  27  10  10  25  27;
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

const subpix3 = [0  7  1  6  3  4  2  5 ;
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

Computes a Peano-Hilbert key for an integer triplet (x,y,z) with x,y,z typically in the range
between 0 and 2^bits-1. Values outside this range can occur when reading across the
periodic box borders.
"""
function peano_hilbert_key(bits::Integer, x::Integer, y::Integer, z::Integer)

    rotation = 1
    key = 0
    mask = 1 << (bits-1)

    @inbounds for _ ∈ 1:bits

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
    get_keylist(h_key::KeyHeader, x0::Array{<:Real}, x1::Array{<:Real})

Get all Peano-Hilbert keys for domain defined by the corner points `x0` and `x1`.
"""
function get_keylist(h_key::KeyHeader, x0::Array{T}, x1::Array{T}) where T

    ix0 = Vector{Int}(undef, 3)
    ix1 = Vector{Int}(undef, 3)
    dix = Vector{Int}(undef, 3)

    nkeys = 1

    @inbounds for i = 1:3
        ix0[i] = get_int_pos( x0[i], h_key.domain_corners[i], h_key.domain_fac, h_key.bits )
        ix1[i] = get_int_pos( x1[i], h_key.domain_corners[i], h_key.domain_fac, h_key.bits )
        dix[i] = ix1[i] - ix0[i] + 1
        nkeys *= dix[i]
    end

    keylist = Vector{Int}(undef, nkeys)

    i = 1
    @inbounds for ix = ix0[1]:ix1[1], iy = ix0[2]:ix1[2], iz = ix0[3]:ix1[3]
        keylist[i] = peano_hilbert_key(h_key.bits, ix, iy, iz)
        i += 1
    end

    return sort!(keylist)
end


"""
    read_positions_from_PH_keys(filebase::String,
                                corner_lowerleft::Array{<:Real}, 
                                corner_upperright::Array{<:Real};
                                parttype::Integer, verbose::Bool)

Finds the Peano-Hilbert keys in a cube defined by `corner_lowerleft` and `corner_upperright`.
"""
function read_positions_from_PH_keys(filebase::String,
                                    corner_lowerleft::Array{<:Real}, 
                                    corner_upperright::Array{<:Real};
                                    parttype::Integer, verbose::Bool)

    # read info blocks once here
    key_info  = read_info(filebase * ".0.key")

    # check if key files are present
    file_key = filebase * ".0.key"
    if !isfile(file_key)
        error("No .key file present! For brute-force read-in set `use_keys=false`")
        flush(stdout)
        flush(stderr)
    end

    if verbose
        @info ".key files found!"
        @info "Calculating peano-hilbert keys..."
        flush(stdout)
        flush(stderr)
        t1 = Dates.now()
    end

    # first read the header
    h_key = read_keyheader(file_key)

    # get a list of the required peano-hilbert keys
    keylist = get_keylist(h_key, corner_lowerleft, corner_upperright)

    if verbose
        t2 = Dates.now()
        @info "$(size(keylist,1)) Peano-Hilbert keys found. Took: $(t2 - t1)"
        @info "Looking for relevant files..."
        flush(stdout)
        flush(stderr)
        t1 = Dates.now()
    end

    # find relevant files
    files = find_files_for_keys(filebase, keylist)
    
    N_files = size(files,1)
    if verbose
        t2 = Dates.now()
        @info "$N_files files found. Took: $(t2 - t1)"

        @info "Searching read positions..."
        flush(stdout)
        flush(stderr)
        println()
        t1 = Dates.now()
    end

    # find all the positions where to read data
    return read_positions_from_keys_files( files, filebase, keylist, key_info; 
                                            parttype, verbose)
end
