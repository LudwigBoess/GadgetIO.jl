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
    global_idxs_to_halo_id(sub_base::String, idxs::Vector{<:Integer}; 
                            parttype::Integer=0)

Converts global halo indices to `HaloID`s.
"""
function global_idxs_to_halo_id(sub_base::String, idxs::AbstractVector{<:Integer}; 
                                parttype::Integer=0)

    # idxs are 0-indexed
    idx_local = idxs .+ 1

    # only works for sorted arrays
    sorted = sortperm(idx_local)
    sort!(idx_local)

    # allocate vector for HaloIDs
    halo_ids = Vector{HaloID}(undef, length(idx_local))

    # if subfind output is only one file you can write the indices directly
    if isfile(sub_base)
        for (i, idx_local) ∈ enumerate(idx_local)
            halo_ids[i] = HaloID(0, idx_local)
        end

    # if there are multiple subfiles we need to arrange the indices
    else
        # get number of sub-files from header
        num_files = read_header(sub_base).num_files

        halos_read = 0
        i = 1
        for filenum = 0:num_files-1

            # read number of halos for this sub-file
            h = read_header(select_file(sub_base, filenum))

            # check currently relevant entries
            sel = findall( halos_read .< idx_local .<= halos_read + h.npart[parttype+1] )

            # write into HaloIDs
            for entry ∈ idx_local[sel]
                # store HaloID
                halo_ids[i] = HaloID(filenum, entry-halos_read)
                # count up entries
                i += 1
            end

            # count up number of halos already read
            halos_read += h.npart[parttype+1]

        end

    end

    halo_ids[sorted]
end


"""
    global_idxs_to_halo_id(sub_base::String, idxs::Vector{<:Integer}; 
                            parttype::Integer=0)

Converts a given number of indices defined by `offset` and `n_to_read` to `HaloID`s.
"""
function global_idxs_to_halo_id(sub_base::String, offset::Integer, n_to_read::Integer; 
                                parttype::Integer=0)

    # read the header 
    h = read_header(sub_base)

    if n_to_read == -1
        finish_read = get_total_particles(h, parttype) - 1
    else
        finish_read = offset + n_to_read
    end

    idxs = offset:finish_read

    return global_idxs_to_halo_id(sub_base, idxs; parttype)
end



"""
    halo_ids_to_read_positions(halo_ids::Vector{HaloID})

Convert a `Vector` of `HaloID`s to a dictionary of `read_positions`. 
To be used with [read_block_filtered](@ref).
"""
function halo_ids_to_read_positions(halo_ids::Vector{HaloID})

    # get all relevant files
    files = unique(getfield.(halo_ids, :file))

    # allocate dict to store IDs per file
    store_arrays = Dict()

    # loop over all halo ids
    for file ∈ files
        
        # filter all HaloIDs in the current file
        sel = findall((x->x.file).(halo_ids) .== file)

        # save all indices
        store_arrays[file] = (p->p.id).(halo_ids[sel])
    end

    # allocate read_positions dict
    read_positions = Dict()

    # loop over file entries
    for file ∈ keys(store_arrays)
        # reduce neighboring block postions
        index, n_to_read = reduce_read_positions(store_arrays[file])

        # store Dicts
        read_positions[file] = Dict( "index" => index, "n_to_read" => n_to_read)
    end

    # finally store all halos to be read
    read_positions["N_part"] = length(halo_ids)

    return read_positions
end


"""
    read_positions_to_halo_ids(read_positions)

Converts `read_positions` to a Vector of [HaloID](@ref)s.
"""
function read_positions_to_halo_ids(read_positions)

    A = Vector{HaloID}(undef, read_positions["N_part"])

    N_read = 1
    for file ∈ keys(read_positions)

        # skip N_part entry
        if file == "N_part"
            continue
        end

        # convert indices and numbers to read to HaloIDs
        for i = 1:length(read_positions[file]["index"]),
            j = 1:read_positions[file]["n_to_read"][i]
            
            # store halo ids
            A[N_read] = HaloID(file, read_positions[file]["index"][i] + j)
            N_read += 1
        end
    end

    return A
end

"""
    IO
"""


"""
    save_halo_ids(filename::String, halo_ids::Vector{HaloID})

Writes a `Vector` of [HaloID](@ref)s to a files.
"""
function save_halo_ids(filename::String, halo_ids::Vector{HaloID})

    # convert halo ids to read positions
    read_positions = halo_ids_to_read_positions(halo_ids)

    # save the read positions
    save_read_positions(filename, read_positions)
end


"""
    save_halo_ids(filename::String)

Loads a `Vector` of [HaloID](@ref)s from a file.
"""
function load_halo_ids(filename::String)
    # load the read positions
    read_positions = load_read_positions(filename)

    # convert read positions to halo ids
    read_positions_to_halo_ids(read_positions)
end
