using ProgressMeter

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
    rvir_block = "RVIR"
    if size(info[getfield.(info, :block_name) .== mass_block],1) == 0
        mass_block = "MTOP"
        rvir_block = "RTOP"
    end

    @showprogress "Reading files..." for i = 0:nfiles-1

        sub_input = select_file(filebase, i)

        M = read_subfind(sub_input, mass_block)
        max_test = findmax(M)

        if max_test[1] > Mmax

            max_file = i
            max_id   = max_test[2][1]

            # store position and virial radius of most massive halo
            POS  = read_subfind(sub_input, "GPOS")[:,max_id]
            RVIR = read_subfind(sub_input, rvir_block)[max_id]

            # store new maximum mass
            Mmax = max_test[1]
        end

    end # for

    return POS, RVIR, HaloID(max_file, max_id)
end




"""
    filter_subfind(sub_base::String, filter_function::Function, files=nothing)

Filters all entries in a subfind file that fulfill the 'filter_funcion' requirements and
returns a `Vector` of [HaloID](@ref)s.

# Examples
```julia
# load packages
using GadgetIO, GadgetUnits

# define filter function
function find_mvir_gt_1e15(filename) 
    h = read_header(filename)
    GU = GadgetPhysical(h) # unit conversion

    # read Mvir and convert to solar masses 
    M = read_subfind(filename, "MVIR") .* GU.m_msun

    return findall(M .> 1.0e15)
end

# basename of subfind output (without .*)
sub_base = /path/to/groups_000/sub_000

# get relevant halos from first 10 files
halo_ids = filter_subfind(sub_base, find_mvir_gt_1e15, 0:9)
```

"""
function filter_subfind(sub_base::String, filter_function::Function, files=nothing)

    # if file range is not defined
    # run over all files
    if isnothing(files)
        h = read_header(sub_base)
        files = 0:h.num_files-1
    end

    # count total number of files
    Nfiles = length(files)

    # allocate storage arrays
    storage_arr = Vector{Vector{Int64}}(undef, Nfiles)

    # loop over all files in parallel
    @showprogress "Filtering files... " for i = 1:Nfiles
        # find correct input file
        sub_input = select_file(sub_base, files[i])

        # apply filter function
        storage_arr[i] = filter_function(sub_input)
    end # for loop

    # count total number of matching entries 
    Nentries = 0
    for i = 1:Nfiles
        Nentries += length(storage_arr[i])
    end

    # allocate empty array for HaloID structs
    A = Vector{HaloID}(undef, Nentries)

    # write HaloIDs to array 
    Ncount = 1
    # loop over all files
    for i = 1:Nfiles
        # loop over all indices where filter function matches
        for j ∈ eachindex(storage_arr[i])
            A[Ncount] = HaloID(files[i], storage_arr[i][j])
            Ncount += 1
        end
    end

    return A
end

