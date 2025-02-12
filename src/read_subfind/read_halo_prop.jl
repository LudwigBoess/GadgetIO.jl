"""
    read_halo_prop(sub_base, blockname::AbstractString, haloid::HaloID; verbose::Bool=true)

Get halo property from block `blockname` by halo id. Returns an array or scalar depending on the block type.
"""
function read_halo_prop(sub_base, blockname::AbstractString, haloid::HaloID; verbose::Bool=true)

    if verbose
        @info "Reading property $blockname of halo $haloid"
        flush(stdout)
        flush(stderr)
    end

    # get full specified block for all halos in file
    sub_input = select_file(sub_base, haloid.file)

    prop = read_subfind(sub_input, blockname, offset=haloid.id - 1, n_to_read=1)

    # prop is either a length=1 vector or a dim×1 matrix, so it is extracted to scalar or vector
    if ndims(prop) == 1
        return prop[1]
    else
        return vec(prop)
    end
end

"""
    read_halo_prop(sub_base, blockname::AbstractString, i_global::Integer; verbose::Bool=true)

Get halo property from block `blockname` by global halo index `i_global` (zero-based index).
"""
function read_halo_prop(sub_base, blockname::AbstractString, i_global::Integer; verbose::Bool=true)
    val, id = read_halo_prop_and_id(sub_base, blockname, i_global; verbose)
    return val
end


"""
    read_halo_prop_and_id(sub_base, blockname::AbstractString, i_global::Integer; verbose::Bool=true)

Get halo property and [`HaloID`](@ref) from block `blockname` by global halo index `i_global` (zero-based index).
`nfiles` should generally be set to `h.num_files`, obtained from `read_header`.
When nfiles is not passed, it is read automatically from the header.
"""
function read_halo_prop_and_id(sub_base, blockname::AbstractString, i_global::Integer; verbose::Bool=true)
    if verbose
        @info "Reading property $blockname of halo at index $i_global"
        flush(stdout)
        flush(stderr)
    end

    # blocks are type specific so we can use this to make our life easier
    parttype = subfind_block_parttype(sub_base, blockname)

    # convert global halo idx to HaloID
    haloid = global_idxs_to_halo_id(sub_base, [i_global]; parttype)[1]

    return read_halo_prop(sub_base, blockname, haloid; verbose), haloid
end



"""
    read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, haloids::AbstractVector{HaloID}; verbose::Bool=true)

Get halo properties defined by an `Array` of blocks for an `Array` of `HaloID`s. Please note that this only works if all blocks are of the same halo type.
    Returns a dictionary with all requested blocks.
"""
function read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, haloids::AbstractVector{HaloID}; verbose::Bool=true)
    if !issorted(haloids)
        @warn "The Vector of HaloIDs is not sorted for requesting the properties from Subfind, the returned properties are returned as if they were sorted, however."
    end

    # check if all blocks are for the same parttype
    parttype = check_subfind_parttype_for_multiple_blocks(sub_base, blocks)

    # get read positions from the array of HaloIDs 
    read_positions = halo_ids_to_read_positions(haloids)

    # read filtered blocks
    return read_blocks_filtered(sub_base, blocks; read_positions, parttype, verbose)
end

"""
    read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, haloid::HaloID; verbose::Bool=true)

Get halo properties defined by an `Array` of blocks for a `HaloID`. Please note that this only works if all blocks are of the same halo type.
Returns a dictionary with all requested blocks.
"""
function read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, haloid::HaloID; verbose::Bool=true)
    # call default function
    read_halo_prop(sub_base, blocks, [haloid]; verbose)
end

"""
    read_halo_prop(sub_base, block::AbstractString, haloids::AbstractVector{HaloID}; verbose::Bool=true)

Get halo properties defined by the block for an `Array` of `HaloID`s.
Returns an array with the requested block.
"""
function read_halo_prop(sub_base, block::AbstractString, haloids::AbstractVector{HaloID}; verbose::Bool=true)
    # call default function
    read_halo_prop(sub_base, [block], haloids; verbose)[block]
end


"""
    read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, i_global::AbstractVector{<:Integer}; verbose::Bool=true)

Get halo properties defined by an `Array` of blocks for an `Array` of global indices. Please note that this only works if all blocks are of the same halo type.
Returns a dictionary with all requested blocks.
"""
function read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, i_global::AbstractVector{<:Integer}; verbose::Bool=true)
    if !issorted(i_global)
        @warn "The Vector of i_global is not sorted for requesting the properties from Subfind, the returned properties are returned as if they were sorted, however."
    end

    # check if all blocks are for the same parttype
    parttype = check_subfind_parttype_for_multiple_blocks(sub_base, blocks)

    # convert global halo idxs to HaloIDs
    haloids = global_idxs_to_halo_id(sub_base, i_global; parttype)

    # get read positions from the array of HaloIDs 
    read_positions = halo_ids_to_read_positions(haloids)

    # read filtered blocks
    return read_blocks_filtered(sub_base, blocks; read_positions, parttype, verbose)
end

"""
    read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, i_global::Integer; verbose::Bool=true)

Get halo properties defined by an `Array` of blocks for a global index. Please note that this only works if all blocks are of the same halo type.
Returns a dictionary with all requested blocks.
"""
function read_halo_prop(sub_base, blocks::AbstractVector{<:AbstractString}, i_global::Integer; verbose::Bool=true)
    # call default function
    read_halo_prop(sub_base, blocks, [i_global]; verbose)
end

"""
    read_halo_prop(sub_base, block::AbstractString, i_global::AbstractVector{<:Integer}; verbose::Bool=true)

Get halo properties defined by the block for an `Array` of global indices.
Returns an array with the requested block.
"""
function read_halo_prop(sub_base, block::AbstractString, i_global::AbstractVector{<:Integer}; verbose::Bool=true)
    # call default function
    read_halo_prop(sub_base, [block], i_global; verbose)[block]
end
