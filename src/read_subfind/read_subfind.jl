"""
    Functions in this file read the subfind output.

"""


using ProgressMeter
using Base.Threads

"""
    struct SubfindHeader

Contains the data of the `HEAD` block in the subfind output

# Fields
| Name                                 | Meaning                                                                                |
| :----------------------------------  | :------------------------------------------------------------------------------------- |
| `nhalos::Int32`                      | number of halos in the output file                                                     |
| `nsubhalos::Int32`                   | number of subhalos in the output file                                                  |
| `nfof::Int32`                        | number of particles in the FoF                                                         |
| `ngroups::Int32`                     | number of large groups in the output file                                              |
| `time::Float64`                      | time / scale factor of the simulation                                                  |
| `z::Float64`                         | redshift of the simulation                                                             |
| `tothalos::UInt32`                   | total number of halos over all output files                                            |
| `totsubhalos::UInt32`                | total number of subhalos over all output files                                         |
| `totfof::UInt32`                     | total number of particles in the FoF                                                   |
| `totgroups::UInt32`                  | 1 if simulation was run with cooling, else 0                                           |
| `num_files::Int32`                   | number of files over which subfind data is distributed                                 |
| `boxsize::Float64`                   | total size of the simulation box                                                       |
| `omega_0::Float64`                   | Omega matter                                                                           |
| `omega_l::Float64`                   | Omega dark enery                                                                       |
| `h0::Float64`                        | little h                                                                               |
| `flag_doubleprecision::Int32`        | 1 if snapshot is in double precision, else 0                                           |
| `flag_ic_info::Int32`                | 1 if initial snapshot file contains an info block, else 0                              |

"""
struct SubfindHeader <: AbstractGadgetHeader
    nhalos::Int32                       # number of halos in the output file
    nsubhalos::Int32                    # number of subhalos in the output file
    nfof::Int32                         # number of particles in the FoF
    ngroups::Int32                      # number of large groups in the output file
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    tothalos::UInt32                    # total number of halos over all output files
    totsubhalos::UInt32                 # total number of subhalos over all output files
    totfof::UInt32                      # total number of particles in the FoF
    totgroups::UInt32                   # total number of large groups over all output files
    num_files::Int32                    # number of files
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32
end

"""
    struct HaloID
        file::Int64
        id::Int64
    end

Stores the subfile that contains the halo in `file` and the position in the block in `id`.
"""
struct HaloID
    file::Int64
    id::Int64
end

"""
    read_subfind_header(filename::String)

Reads the header of a subfind file or file base (without .0, .1, etc.) into a SubfindHeader struct.
"""
function read_subfind_header(filename::String)

    filename = select_file(filename, 0)
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
            println("incorrect file format encountered when reading header of", filename)
        end
    end

    #println("Reading snapshot format: ", snap_format)

    if snap_format == 2
        seek(f, 16)
        skip_line = read(f, Int32)
    end

    nhalos = read(f, Int32)
    nsubhalos = read(f, Int32)
    nfof = read(f, Int32)
    ngroups = read(f, Int32)

    # rest of numbers and massarr are obsolete
    DUMMY = read!(f, Array{Int32,1}(undef,2))
    DUMMY = read!(f, Array{Float64,1}(undef,6))

    time = read(f, Float64)
    z = read(f, Float64)

    # flags not needed
    DUMMY = read(f, Int32)
    DUMMY = read(f, Int32)

    tothalos = read(f, UInt32)
    totsubhalos = read(f, UInt32)
    totfof = read(f, UInt32)
    totgroups = read(f, UInt32)

    DUMMY = read!(f, Array{Int32,1}(undef,2))
    DUMMY = read(f, Int32)

    num_files = read(f, Int32)
    boxsize = read(f, Float64)
    omega_0 = read(f, Float64)
    omega_l = read(f, Float64)
    h0 = read(f, Float64)

    DUMMY = read(f, Int32)
    DUMMY = read(f, Int32)
    DUMMY = read!(f, Array{UInt32,1}(undef,6))
    DUMMY = read(f, Int32)

    flag_doubleprecision = read(f, Int32)
    flag_ic_info = read(f, Int32)

    close(f)

    return SubfindHeader(nhalos, nsubhalos, nfof, ngroups,
                         time, z,
                         tothalos, totsubhalos, totfof, totgroups,
                         num_files,
                         boxsize, omega_0, omega_l, h0,
                         flag_doubleprecision,
                         flag_ic_info)
end

"""
    read_subfind(filename::String, blockname::String; return_haloid::Bool=false)

Reads a block of a subfind file.

If `return_haloid` is `true`, returns a tuple of the block array and the corresponding `HaloID`s.
"""
function read_subfind(filename::String, blockname::String; return_haloid::Bool=false)
    
    # read the info block
    info = read_info(filename)

    blocknames = getfield.(info, :block_name)
    ind = findfirst(==(blockname), blocknames)

    # the the block is not contained in the file throw and error
    if isnothing(ind)
        error("Block $blockname not present!")
    end

    # get the relevant entry
    info_selected = info[ind]

    # blocks are type specific so we can use this to make our life easier
    parttype = findfirst(==(1), info_selected.is_present) - 1

    return read_block(filename, blockname; info=info_selected, parttype, return_haloid)
end


"""
    read_subfind(filename::String, blockname::String, ids::AbstractVector{<:Integer}; return_haloid::Bool=false)

Reads the block at the given subfind ids (0-indexed).

If `return_haloid` is `true`, returns a tuple of the block array and the corresponding `HaloID`s.
"""
function read_subfind(filename::String, blockname::String, ids::AbstractVector{<:Integer}; return_haloid::Bool=false)
    # shift to 1-indexed
    ind_ids = ids .+ 1
    
    # read the info block
    info = read_info(filename)

    blocknames = getfield.(info, :block_name)
    ind = findfirst(==(blockname), blocknames)

    # the the block is not contained in the file throw and error
    if isnothing(ind)
        error("Block $blockname not present!")
    end

    # get the relevant entry
    info_selected = info[ind]

    # blocks are type specific so we can use this to make our life easier
    parttype = findfirst(==(1), info_selected.is_present) - 1


    if isfile(filename)
        if return_haloid
            block, haloids = read_block(filename, blockname; inf=info_selected, parttype, return_haloid)
            if ndims(block) == 1
                return block[ind_ids], haloids
            else
                return block[:,ind_ids], haloids
            end
        else
            block = read_block(filename, blockname; inf=info_selected, parttype)
            if ndims(block) == 1
                return block[ind_ids]
            else
                return block[:,ind_ids]
            end
        end
    end


    h = read_header(filename)
    maxind = maximum(ind_ids)

    # initialize arrays holding blocks and haloids
    block_arr = []
    if return_haloid
        haloids_arr = []
    end

    # loop over all sub-files
    n_already_read_in = 0
    @inbounds for i = 0:(h.num_files-1)
        if return_haloid
            block, haloids = read_block(filename*".$i", blockname; parttype, return_haloid, thisfilenum=i)
            push!(haloids_arr, haloids)
        else
            block = read_block(filename*".$i", blockname; parttype)
        end
        push!(block_arr, block)

        n_already_read_in += length(block)

        if n_already_read_in ≥ maxind
            break
        end
    end

    if return_haloid
        return get_lazy_vcat_index.((block_arr,), ind_ids), get_lazy_vcat_index.((haloids_arr,), ind_ids)
    else
        return get_lazy_vcat_index.((block_arr,), ind_ids)
    end
end

"""
    get_lazy_vcat_index(arr::AbstractVector, ind)

For an array of arrays `[a, b, c]`, where `a`, `b`, and `c` are arrays, returns the value of `vcat(a, b, c)[ind]`.

This is not exported.
"""
function get_lazy_vcat_index(arr::AbstractVector, ind)
    @assert ind ≥ 0

    for a in arr
        n = length(a)
        if ind ≤ n
            return a[ind]
        end

        ind -= n
    end

    throw(BoundsError("attempt to access lazy array of arrays at too large index"))
end
