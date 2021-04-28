using Dates
using Base.Threads

"""
    find_read_positions(files::Array{<:Integer}, filebase::String, 
                             blocks::Array{String}, parttype::Integer)

Helper function to get positions and length of particle blocks in files.
"""
function find_read_positions(files::Array{<:Integer}, filebase::String, 
                             blocks::Array{String}, parttype::Integer,
                             keylist::Array{<:Integer}, key_info::Array{InfoLine},
                             verbose::Bool)

    # store number of file
    N_files = size(files,1)

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
        keys_in_file = read_block(filename_keyfile, "KEY",
                                          info = key_info[getfield.(key_info, :block_name) .== "KEY"][1],
                                          parttype = parttype)

        if verbose
            @info "Calculating index list..."
            t1 = Dates.now()
        end

        index_list = get_index_list(keylist, keys_in_file)

        if verbose
            t2 = Dates.now()
            @info "Index list done. Took: $(t2 - t1)"
            @info "Reading $(size(index_list,1)) key segments..."
        end

        # number of particles associated with PH key
        part_per_key = read_block(filename_keyfile, "NKEY",
                                          info = key_info[getfield.(key_info, :block_name) .== "NKEY"][1],
                                          parttype = parttype)


        # offsets in the blocks to get to the relevant particles
        offset_key = read_block(filename_keyfile, "OKEY",
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
            @info "Reduced independent blocks from $(size(offset_key,1)) to $(size(use_block[use_block],1))"
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

    end # for i = 1:size(files,1)

    return file_offset_key, file_part_per_key, file_block_positions
end

function get_index_bounds(ids::Vector{<:Integer}, low_bounds::Vector{<:Integer}, high_bounds::Vector{<:Integer})

    nids = size(ids,1)
    nbounds = size(low_bounds,1)

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

                    if icountids > nids
                        lend2 = true
                    end # icountids > nids
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

                            if icountids > nids
                                lend2 = true
                            end
                        end # if ids[icountids] >= low_bounds[icountbounds]

                    end # while !lend2

                end # if icountids <= nids

            elseif ids[icountids] > high_bounds[icountbounds] # if ids[icountids] < low_bounds[icountbounds]
                icountbounds += 1

                if icountbounds <= nbounds
                    lend2 = false

                    while !lend2
                        if ids[icountids] <= high_bounds[icountbounds]
                            lend2 = true
                        else # ids[icountids] <= high_bounds[icountbounds]
                            icountbounds += 1
                            if icountbounds > nbounds
                                lend2 = true
                            end # icountbounds > nbounds
                        end # ids[icountids] <= high_bounds[icountbounds]
                    end # while !lend2

                end # icountbounds <= nbounds

            end # if ids[icountids] < low_bounds[icountbounds]

        end # if low_bounds[icountbounds] <= ids[icountids] <= high_bounds[icountbounds]

        if icountids > nids
            lend = true
        end
        if icountbounds > nbounds
            lend = true
        end

        #
        # println("icountids     = ", icountids)
        # println("icountbounds  = ", icountbounds)

    end # while !lend

    if icountall > 1
        ind_out = ind_all[1:icountall-1]
        return ind_out
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

    mask = falses(size(low_list,1))

    for key in keylist
        @. mask = mask | ( (key >= low_list ) & ( key <= high_list ))
    end

    return Int64.(unique!(sort!(file_list[mask])))
end


"""
    get_index_list(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

Finds the indices at which `list_to_check` contains elements from `list_to_find`.
If both either of the lists are not sorted it uses a `Dict` lookup, otherwise it uses a `Array` forward-search.
"""
function get_index_list(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

    # check if the lists are sorted
    if ( issorted(list_to_find) && issorted(list_to_check) )
        return get_index_list_arr(list_to_find, list_to_check)
    else
        return get_index_list_set(list_to_find, list_to_check)
    end
end

"""
    get_index_list_arr(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`. Uses forward-searching in sorted array.
"""
@inline function get_index_list_arr(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

    narr1 = size(list_to_find,1)
    narr2 = size(list_to_check,1)

    ind_all   = zeros(Int64, narr1)

    icountall     = 1
    icountarr1    = 1
    icountarr2    = 1

    lend = false

    while !lend

        if list_to_check[icountarr2] == list_to_find[icountarr1]

            ind_all[icountall] = icountarr2
            icountall  += 1
            icountarr1 += 1
            icountarr2 += 1
        else  # list_to_check[icountnotarr2] == list_to_find[icountnotar1]
            if list_to_check[icountarr2] < list_to_find[icountarr1]
                icountarr2    += 1
            else # list_to_check[icountnotarr2] < list_to_find[icountnotar1]
                icountarr1    += 1
            end # list_to_check[icountnotarr2] < list_to_find[icountnotar1]
        end # list_to_check[icountnotarr2] == list_to_find[icountnotar1]
        if (icountarr2 > narr2 ) || (icountarr1 > narr1)
            lend = true
        end
    end # while !lend

    if icountall > 1
        return ind_all[1:icountall-1]
    end
end

"""
    get_index_list_arr(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`. Uses forward-searching in sorted array.
"""
@inline function get_index_list_arr_sort(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

    narr1 = size(list_to_find,1)
    narr2 = size(list_to_check,1)

    ind_all   = zeros(Int64, narr1)

    icountall     = 1
    icountarr1    = 1
    icountarr2    = 1

    iiarr2 = sortperm(idarr2[:,1], alg=QuickSort)

    lend = false

    while !lend

        if list_to_check[iiarr2[icountarr2]] == list_to_find[icountarr1]

            ind_all[icountall] = iiarr2[icountarr2]
            icountall  += 1
            icountarr1 += 1
            icountarr2 += 1
        else  # list_to_check[icountnotarr2] == list_to_find[icountnotar1]
            if list_to_check[iiarr2[icountarr2]] < list_to_find[icountarr1]
                icountarr2    += 1
            else # list_to_check[icountnotarr2] < list_to_find[icountnotar1]
                icountarr1    += 1
            end # list_to_check[icountnotarr2] < list_to_find[icountnotar1]
        end # list_to_check[icountnotarr2] == list_to_find[icountnotar1]
        if (icountarr2 > narr2 ) || (icountarr1 > narr1)
            lend = true
        end
    end # while !lend

    if icountall > 1
        return ind_all[1:icountall-1]
    end
end

"""
    get_index_list_dict(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`. Uses a `Dict` for lookup -> slower than the array search, but works on unsorted arrays.
"""
@inline function get_index_list_dict(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

    dict = Dict((n, i) for (i, n) in enumerate(list_to_check))
    result = Vector{Int}(undef, size(list_to_find,1))
    len = 0

    for k in list_to_find
        i = get(dict, k, nothing)
        if i !== nothing
            len += 1
            @inbounds result[len] = i
        end
    end
    return resize!(result, len)
end

"""
    get_index_list_set(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`.
Uses a `Set` for lookup -> slower than the array search, but works on unsorted arrays
like `Dict`, just faster.
"""
@inline function get_index_list_set(list_to_find::Array{<:Integer}, list_to_check::Array{<:Integer})

    set = Set(list_to_check)

    return findall(in.(list_to_find, (set,)))
end


"""
   function join_blocks(offset_key, part_per_key)
    
Joins neigboring blocks to simplify read-in.
"""
@inline function join_blocks(offset_key, part_per_key)

    use_block = trues(size(offset_key,1))

    icount = 1

    for i = 2:size(offset_key,1)-1

        if offset_key[i] == offset_key[icount] + part_per_key[icount]
            part_per_key[icount] += part_per_key[i]
            use_block[i] = false
        else
            icount += 1
        end # if

    end # for

    return use_block, part_per_key
end
