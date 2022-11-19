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
function global_idxs_to_halo_id(sub_base::String, idxs::Vector{<:Integer}; 
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

    idxs = collect(offset:finish_read)

    return global_idxs_to_halo_id(sub_base, idxs; parttype)
end