"""
Batched multi-block subfind reader for the merger-tree pipeline.

`load_subhalo_catalogue` and `load_subfind_meta` each call `read_subfind`
many times — once per block. For a multi-file subfind output (e.g. 16
sub-files), every call re-opens every sub-file, re-reads its header,
re-reads INFO, and re-scans block positions. On NFS/Lustre this is
IO-bound and dominates runtime.

`read_subfind_blocks(sub_base, blocks; parttype)` opens each sub-file
once, caches the header/info/block-position lookup, and reads every
requested block in a single pass. Threaded over sub-files.

Inputs: `blocks` must all belong to the same `parttype`. The result is
a `Dict{String,Array}` keyed by block name.
"""

# ---------------------------------------------------------------------------
"""
    read_subfind_blocks(sub_base, blocks; parttype) -> Dict{String,Array}

Read multiple subfind blocks in one batched pass over the sub-files.
All requested blocks must share the same `parttype` (e.g. all subhalo
blocks or all FOF-halo blocks).

The result vectors/matrices are concatenated across sub-files in sub-file
order, identical to what `read_subfind(sub_base, block)` would return for
each block individually.
"""
function read_subfind_blocks(sub_base::AbstractString,
                             blocks::AbstractVector{<:AbstractString};
                             parttype::Integer)

    isempty(blocks) && return Dict{String,Any}()

    # discover sub-files
    h_outer = read_header(sub_base)
    num_files = isfile(sub_base) ? 1 : Int(h_outer.num_files)

    # per-file headers + npart of the requested parttype
    headers = Vector{SnapshotHeader}(undef, num_files)
    nparts  = Vector{Int64}(undef, num_files)
    for fi in 1:num_files
        fn = select_file(sub_base, fi - 1)
        headers[fi] = read_header(fn)
        nparts[fi]  = headers[fi].npart[parttype + 1]
    end
    n_total = sum(nparts)

    # info from the first file that actually contains this parttype
    info_list = nothing
    for fi in 1:num_files
        if nparts[fi] > 0
            info_list = read_info(select_file(sub_base, fi - 1))
            break
        end
    end
    info_list === nothing && error("No sub-file contains parttype $parttype")

    # per-block InfoLine (same across sub-files in a well-formed dataset)
    block_info = Dict{String,InfoLine}()
    for b in blocks
        bname = String(strip(b))
        idx = findfirst(info_list) do il
            String(strip(il.block_name)) == bname
        end
        idx === nothing && error("Block $bname not present in subfind output")
        block_info[bname] = info_list[idx]
    end

    # allocate output arrays
    out = Dict{String,Any}()
    for b in blocks
        bname = String(strip(b))
        info = block_info[bname]
        if info.n_dim == 1
            out[bname] = Vector{info.data_type}(undef, n_total)
        else
            out[bname] = Matrix{info.data_type}(undef, info.n_dim, n_total)
        end
    end

    # cumulative offsets (start of each sub-file in the global array)
    offsets = Vector{Int64}(undef, num_files)
    acc = Int64(0)
    for fi in 1:num_files
        offsets[fi] = acc
        acc += nparts[fi]
    end

    # threaded read: one sub-file per task, all blocks at once
    Threads.@threads for fi in 1:num_files
        nparts[fi] == 0 && continue

        fn = select_file(sub_base, fi - 1)
        h_local = headers[fi]
        positions = get_block_positions(fn)
        n_local = Int(nparts[fi])
        my_offset = offsets[fi]

        io = open(fn, "r")
        try
            for b in blocks
                bname = String(strip(b))
                info = block_info[bname]
                arr = out[bname]
                read_block!(arr, io, 0, my_offset, n_local;
                            parttype,
                            block_position = positions[bname],
                            info, h = h_local)
            end
        finally
            close(io)
        end
    end

    return out
end
