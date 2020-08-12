"""
    Functions in this file read the subfind output.

"""

struct SubfindHeader
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
    num_colors::Int32                   # number of colors
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32
end

function read_subfind_header(filename::String)

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


function read_subfind(filename::String, blockname::String)

    info = read_info(filename)

    info = info[getfield.(info, :block_name) .== blockname][1]

    parttype = findall(info.is_present .== 1)[1] - 1

    return read_snap(filename, blockname, parttype)
end


using ProgressMeter
using Base.Threads

struct SubfindID
    file::Int64
    id::Int64
end

@inline function select_file(filebase::String, i::Int, nfiles::Int)

    if nfiles > 1
        sub_input = filebase * ".$i"
        if !isfile(sub_input)
            error("File $sub_input not present!")
        end
    else  # nfiles > 1
        sub_input = filebase

        if !isfile(sub_input)
            sub_input = filebase * ".0"
        end

        if !isfile(sub_input)
            error("Subfind file not present!")
        end

    end  # nfiles > 1

    return sub_input
end


"""
    find_most_massive_halo(filebase::String [, nfiles::Int=1])

Reads the selected file and its subfiles and returns position, virial radius
and a SubfindID object that contains the subfile which contains the most massive
halo and the position in the block.

"""
function find_most_massive_halo(filebase::String, nfiles::Int=1)

    Mmax     = 0.0
    POS      = zeros(Float32, 3)
    RVIR     = 0.0
    max_file = 0
    max_id   = 0

    @showprogress "Reading files..." for i = 0:nfiles-1

        sub_input = select_file(filebase, i, nfiles)

        M = read_subfind(sub_input, "MVIR")
        max_test = findmax(M)

        if max_test[1] > Mmax

            max_file = i
            max_id   = max_test[2][1]

            # store position and virial radius of most massive halo
            POS  = read_subfind(sub_input, "GPOS")[max_id,:]
            RVIR = read_subfind(sub_input, "RVIR")[max_id]

            # store new maximum mass
            Mmax = max_test[1]
        end

    end # for

    return POS, RVIR, SubfindID(max_file, max_id)
end

struct SubfindFilter
    file::Int64
    index::Array{Int64}
end


"""
    filter_subfind(filebase::String, blockname::String, filter_function::Function [, nfiles::Integer=1])

Selects entries in subfind block that fulfill the 'filter_funcion' requirements and
returns a 'SubfindFilter' object.

# Examples
```jldoctest
julia> find_mass_gt_1e15(M) = ( (M > 1.e15) ? true : false )
find_mass_gt_1e15 (generic function with 1 method)
julia> filtered_subfind = filter_subfind(filebase, "MVIR", find_mass_gt_1e15)
[...]
```

"""
function filter_subfind(filebase::String, blockname::String, filter_function::Function, nfiles::Integer=1)

    # allocate empty array for SubfindFilter objects
    A = Vector{SubfindFilter}(undef, 0)

    # loop over all files in parallel
    @threads for i = 0:nfiles-1

        sub_input = select_file(filebase, i, nfiles)

        block = read_subfind(sub_input, blockname)

        selection = filter_function.(block)
        correct_selection = findall(selection)

        if length(selection[selection]) > 0

            # create array of integers for easy storing
            id_array = collect(1:length(selection))[correct_selection]

            # create SubfindFilter object and push it to the storage array
            push!(A, SubfindFilter(i, id_array))

        end # if entries > 0

    end # for loop

    return A
end

function read_particles_in_halo(sub_file::String, snap_filebase::String,
                                halo_id::Int64, blocks::Vector{String};
                                parttype::Integer=0, Rhalf_factor::AbstractFloat=5.0)


    error("Not finished!")

    # read subfind data

    pos      = read_subfind(sub_input, "SPOS")[halo_id,:]
    rhalf    = read_subfind(sub_input, "RHMS")[halo_id]

    sub_info = read_info(sub_file)
    pos = get_block_positions(sub_file)
    ids = read_block_with_offset(sub_info, "KEY",
                                      info = sub_info[getfield.(sub_info, :block_name) .== "PID"][1],
                                      parttype = 0)
end
