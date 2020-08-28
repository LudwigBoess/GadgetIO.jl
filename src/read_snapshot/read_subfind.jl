"""
    Functions in this file read the subfind output.

"""


using ProgressMeter
using Base.Threads

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

Reads the header of a subfind file into a SubfindHeader struct.
"""
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

"""
    read_subfind(filename::String, blockname::String)

Reads a block of a subfind file.
"""
function read_subfind(filename::String, blockname::String)

    # read the info block
    info = read_info(filename)

    # the the block is not contained in the file throw and error
    if length(info[getfield.(info, :block_name) .== blockname]) == 0
        error("Block $blockname not present!")
    end

    # get the relevant entry
    info_selected = info[getfield.(info, :block_name) .== blockname][1]

    # blocks are type specific so we can use this to make our life easier
    parttype = findall(info_selected.is_present .== 1)[1] - 1

    return read_block_by_name(filename, blockname, info = info_selected, parttype = parttype)
end


"""
    find_most_massive_halo(filebase::String [, nfiles::Int=1])

Reads the selected file and its subfiles and returns position, virial radius
and a HaloID object that contains the subfile which contains the most massive
halo and the position in the block.

"""
function find_most_massive_halo(filebase::String, nfiles::Int=1)

    Mmax     = 0.0
    POS      = zeros(Float32, 3)
    RVIR     = 0.0
    max_file = 0
    max_id   = 0

    sub_input = select_file(filebase, 0)
    # read the info block
    info = read_info(sub_input)

    # check if the massblock is called MVIR or MTOP
    mass_block = "MVIR"
    if length(info[getfield.(info, :block_name) .== mass_block]) == 0
        mass_block = "MTOP"
    end

    @showprogress "Reading files..." for i = 0:nfiles-1

        sub_input = select_file(filebase, i)

        M = read_subfind(sub_input, mass_block)
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

    return POS, RVIR, HaloID(max_file, max_id)
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
    A = Vector{HaloID}(undef, 0)

    # loop over all files in parallel
    @threads for i = 0:nfiles-1

        # find correct input file
        sub_input = select_file(filebase, i)

        # read block
        block = read_subfind(sub_input, blockname)

        # apply filter function
        selection = filter_function.(block)

        # select matches
        correct_selection = findall(selection)

        if length(selection[selection]) > 0

            # create array of integers for easy storing
            id_array = collect(1:length(selection))[correct_selection]

            # create HaloID object and push it to the storage array
            for j = 1:length(id_array)
                push!(A, HaloID(i, id_array[j]))
            end

        end # if entries > 0

    end # for loop

    return A
end

