using Dates
using Base.Threads

"""
    read_positions_from_PH_keys(files::Vector{<:Integer}, filebase::String, 
                                blocks::AbstractVector{String}, parttype::Integer,
                                keylist::Vector{<:Integer}, key_info::Vector{InfoLine},
                                verbose::Bool)

Helper function to get positions and length of particle blocks in files.
"""
function read_positions_from_keys_files(files::Vector{<:Integer}, filebase::String,
                                        keylist::Vector{<:Integer}, key_info::Vector{InfoLine};
                                        parttype::Integer, verbose::Bool)

    # store number of file
    N_files = length(files)

    # allocate arrys to store reading information
    d = Dict()

    d["N_part"] = 0

    @inbounds for i = 1:N_files

        filename = select_file(filebase, files[i])

        # read header of the file
        h = read_header(filename)

        # read header of the keyfile
        filename_keyfile = filename * ".key"
        key_h     = read_header(filename_keyfile)

        # overwrite number of files to re-use simple block reading
        key_h.num_files = 1

        if iszero(h.npart[parttype+1])
            @info "No particles of type $parttype in file $(i)!"
            continue
        end

        # Vector of field names (KEY, NKEY, OKEY)
        fields = getfield.(key_info, :block_name)

        # Reads the following in, respectively:
        # - key file data
        # - number of particles associated with PH key
        # - offsets in the blocks to get to the relevant particles
        keys_in_file = read_block(filename_keyfile, "KEY", 
                                   info = key_info[findfirst(==("KEY"), fields)],
                                   h = key_h; 
                                   parttype)

        part_per_key = read_block(filename_keyfile, "NKEY", 
                                   info = key_info[findfirst(==("NKEY"), fields)],
                                   h = key_h;
                                   parttype)


        offset_key   = read_block(filename_keyfile, "OKEY", 
                                   info = key_info[findfirst(==("OKEY"), fields)],
                                   h = key_h; 
                                   parttype)


        if verbose
            @info "Calculating index list..."
            t1 = Dates.now()
        end

        index_list = get_index_list(keylist, keys_in_file)

        if verbose
            t2 = Dates.now()
            @info "Index list done. Took: $(t2 - t1)"
            @info "Reading $(length(index_list)) key segments..."
        end

        # sort the offset arrays for simplification step
        sorted_offset = sortperm(offset_key[index_list])

        # save sorted arrays
        offset_key   = @views offset_key[index_list[sorted_offset]]
        part_per_key = @views part_per_key[index_list[sorted_offset]]

        # check if blocks can be joined
        use_block, part_per_key = join_blocks(offset_key, part_per_key)

        if verbose
            @info "Reduced independent blocks from $(length(offset_key)) to $(count(use_block))"
        end

        # store the arrays for later reading
        d[i]   = Dict("index"     => offset_key[use_block],
                      "n_to_read" => part_per_key[use_block])

        # sum up all particles to read
        d["N_part"] += sum(d[i]["n_to_read"])

    end # for i = 1:N_files

    return d
end

"""
    read_key_block( filename_keyfile::String, key_info::Vector{InfoLine}, 
                    fields::Vector{String}, key::String, parttype::Integer)

Returns key block for different `key`. Valid values of `key`: KEY, NKEY, OKEY
"""
function read_key_block(filename_keyfile::String, key_info::Vector{InfoLine}, 
                        fields::Vector{String}, key::String, parttype::Integer)
    read_block(filename_keyfile, key,
               info = key_info[findfirst(==(key), fields)],
               parttype = parttype)
end


"""
    get_index_bounds(ids::Vector{<:Integer}, low_bounds::Vector{<:Integer}, high_bounds::Vector{<:Integer})

Returns sorted `Vector` of indices `i`, for which `low_bounds[i] ≤ ids[j] ≤ high_bounds[i]` for any `j`.
All parameters `ids`, `low_bounds`, and `high_bounds` have to already be sorted.
"""
function get_index_bounds(ids::Vector{<:Integer}, low_bounds::Vector{<:Integer}, high_bounds::Vector{<:Integer})
    
    nids    = length(ids)
    nbounds = length(low_bounds)
    ind_all = Vector{Int}(undef, min(nids, nbounds))

    icountall = 1
    icountids = 1
    @inbounds for icountbounds = 1:length(low_bounds)

        while icountids ≤ nids && ids[icountids] < low_bounds[icountbounds]
            icountids += 1
        end

        if icountids ≤ nids && ids[icountids] ≤ high_bounds[icountbounds]
            ind_all[icountall] = icountbounds
            icountall += 1
            icountids += 1
        end

    
    end

    return ind_all[1:(icountall-1)]
end





"""
    find_files_for_keys(filebase::String, nfiles::Integer, keylist::Vector{<:Integer})

Selects the files in which the particles associated with the given Peano-Hilbert keys are stored.
"""
function find_files_for_keys(filebase::String, keylist::Vector{<:Integer})

    file_key_index = filebase * ".key.index"

    # if index file does not exist all key files need to be read
    if !isfile(file_key_index)
        return error("No .key.index file present!")
    end

    # get the data from the index file (low and high lists are sorted)
    low_list, high_list, file_list = read_key_index(file_key_index)

    # get files that contain the wanted particles (parameters have to be sorted)
    index_bounds = get_index_bounds(sort(keylist), low_list, high_list)

    file_sort = file_list[index_bounds] |> unique! |> sort!

    return Int64.(file_sort)
end


"""
    get_index_list(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

Finds the indices at which `list_to_check` contains elements from `list_to_find`.
If both either of the lists are not sorted it uses a `Set` lookup, otherwise it uses a `Vector` forward-search.
"""
function get_index_list(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

    # check if the lists are sorted
    if ( issorted(list_to_find) && issorted(list_to_check) )
        return get_index_list_arr(list_to_find, list_to_check)
    else
        return get_index_list_dict(list_to_find, list_to_check)
    end
end

"""
    get_index_list_arr(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`. Uses forward-searching in sorted array.
Both arrays have to be sorted.
"""
@inline function get_index_list_arr(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})
    if isempty(list_to_find) || isempty(list_to_check)
        return Int64[]
    end

    narr1 = length(list_to_find)
    narr2 = length(list_to_check)

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
    else
        return Int64[]
    end
end


"""
    get_index_list_dict(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`. Uses a `Dict` for lookup -> slower than the array search, but works on unsorted arrays.
"""
@inline function get_index_list_dict(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

    dict = Dict((n, i) for (i, n) in enumerate(list_to_check))
    result = Vector{Int}(undef, length(list_to_find))
    len = 0

    for k in list_to_find
        i = get(dict, k, nothing)
        if !isnothing(i)
            len += 1
            @inbounds result[len] = i
        end
    end
    return resize!(result, len)
end

"""
    get_index_list_set(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

Get positions in `list_to_check` where `list_to_check` matches `list_to_find`.
Uses a `Set` for lookup -> slower than the array search `get_index_list_arr`, but works on unsorted arrays
like `get_index_list_dict`, just faster.
"""
@inline function get_index_list_set(list_to_find::Vector{<:Integer}, list_to_check::Vector{<:Integer})

    set = Set(list_to_find)

    return findall(in.(list_to_check, (set,)))
end


"""
   join_blocks(offset_key, part_per_key)
    
Joins neighboring blocks to simplify read-in.
"""
@inline function join_blocks(offset_key::AbstractVector{<:Integer}, part_per_key::AbstractVector{<:Integer})

    use_block = trues(length(offset_key))

    icount = 1

    @inbounds for i = 2:length(offset_key)-1

        if offset_key[i] == offset_key[icount] + part_per_key[icount]
            part_per_key[icount] += part_per_key[i]
            use_block[i] = false
        else
            icount += 1
        end # if

    end # for

    return use_block, part_per_key
end
