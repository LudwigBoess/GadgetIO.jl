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
    if size(info[getfield.(info, :block_name) .== mass_block])[1] == 0
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
            POS  = read_subfind(sub_input, "GPOS")[max_id,:]
            RVIR = read_subfind(sub_input, rvir_block)[max_id]

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

        if size(selection[selection])[1] > 0

            # create array of integers for easy storing
            id_array = collect(1:size(selection)[1])[correct_selection]

            # create HaloID object and push it to the storage array
            for j = 1:size(id_array)[1]
                push!(A, HaloID(i, id_array[j]))
            end

        end # if entries > 0

    end # for loop

    return A
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
function filter_subfind(filebase::String, filter_function::Function, nfiles::Integer=1)

    # allocate empty array for SubfindFilter objects
    A = Vector{HaloID}(undef, 0)

    # loop over all files in parallel
    @threads for i = 0:nfiles-1

        # find correct input file
        sub_input = select_file(filebase, i)

        # apply filter function
        id_array = filter_function(sub_input)

        if size(id_array,1) > 0
            # create HaloID object and push it to the storage array
            for j = 1:size(id_array,1)
                push!(A, HaloID(i, id_array[j]))
            end

        end # if entries > 0

    end # for loop

    return A
end


