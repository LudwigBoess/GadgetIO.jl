"""
    read_subfind(filename::String, blockname::String)

# Example
To read e.g. the virial radius of halos use
```julia
R_vir = read_subfind(filename, "RVIR")
```
"""
function read_subfind(filename::String, blockname::String; 
                      info::Union{Nothing,InfoLine}=nothing,
                      offset::Integer=0, n_to_read::Integer=-1,           
                      return_haloid::Bool=false)
    

    # blocks are type specific so we can use this to make our life easier
    parttype = subfind_block_parttype(filename, blockname)

    block = read_block(filename, blockname; parttype, info, offset, n_to_read)

    if return_haloid
        return block, global_idxs_to_halo_id(filename, offset, n_to_read; 
                                             parttype)
    else
        return block
    end
end

"""
    read_subfind(filename::String, blockname::String, ids::AbstractVector{<:Integer}; return_haloid::Bool=false)

Reads the block at the given subfind indices (`ids`, 0-indexed).

If `return_haloid` is `true`, returns a tuple of the block array and the corresponding `HaloID`s.

# Example
```julia
# Read the virial masses of the first four halos in subfind:
mvir = read_subfind(filebase, "MVIR", [0, 1, 2, 3])

# or:
mvir, haloids = read_subfind(filebase, "MVIR", [0, 1, 2, 3]; return_haloid=true)
```
"""
function read_subfind(filename::String, blockname::String, ids::AbstractVector{<:Integer}; 
                        info::Union{Nothing,InfoLine}=nothing,
                        offset::Integer=0, n_to_read::Integer=-1,           
                        return_haloid::Bool=false)
    # shift to 1-indexed
    ind_ids = ids .+ 1
    
    # read the block data
    block = read_subfind(filename, blockname; info, offset, n_to_read)

    if ndims(block) == 1
        block = block[ind_ids]
    else
        block = block[:,ind_ids]
    end

    if return_haloid
        # blocks are type specific so we can use this to make our life easier
        parttype = subfind_block_parttype(filename, blockname)
        return block, global_idxs_to_halo_id(filename, ids; parttype)
    else
        return block
    end

end

"""
    subfind_block_parttype(filename, blockname)

Get the particle type for which a subfind block is relevant.
"""
function subfind_block_parttype(filename, blockname, info=nothing)

    if isnothing(info)
        # read the info block
        info = check_info(filename, blockname)
    end

    # read the info block
    info = check_info(filename, blockname)

    # blocks are type specific so we can use this to make our life easier
    return findfirst(==(1), info.is_present) - 1
end

"""
    check_subfind_parttype_for_multiple_blocks(sub_base, blocks::Vector{AbstractString})

Checks if all requested blocks are for the same halo type and returns the correct halo type if true.
"""
function check_subfind_parttype_for_multiple_blocks(sub_base, blocks::Vector{AbstractString})

    # blocks are type specific so we can use this to make our life easier
    parttype = subfind_block_parttype(sub_base, blocks[1])

    # check if all requested blocks are for the same halo type
    for block ∈ blocks 
        
        # read parttype of current block
        parttype_block = subfind_block_parttype(sub_base, block)

        # error handling
        if parttype_block != parttype 
            error("All requested blocks must be for the same halo type. Block $block is not available for halo type $parttype.")
        end
    end

    # if all parttypes match we can return the correct parttype
    return parttype
end

"""
    read_subfind_length(filename::String, blockname::String)

Get number of entries for block `blockname`. Uses the header, so it is faster than reading the whole block.
"""
function read_subfind_length(filename::String, blockname::String)
    
    # blocks are type specific so we can use this to make our life easier
    parttype = subfind_block_parttype(filename, blockname)

    # read header
    h = head_to_struct(filename)

    return h.npart[parttype+1]
end
