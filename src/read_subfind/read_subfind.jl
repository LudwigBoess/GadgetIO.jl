"""
    read_subfind(filename::String, blockname::String)

Reads a block of a subfind file.

If `return_haloid` is `true`, returns a tuple of the block array and the corresponding `HaloID`s.
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
    read_subfind_length(filename::String, blockname::String)

Get number of entries for block `blockname`. Uses the header, so it is faster than reading the whole block.
"""
function read_subfind_length(filename::String, blockname::String)
    
    # blocks are type specific so we can use this to make our life easier
    parttype = subfind_block_parttype(filename, blockname)

    h = head_to_obj(filename)

    return h.npart[parttype+1]
end



"""
    get_lazy_vcat_indices(arr::AbstractVector, inds)

For an array of arrays `[a, b, c]`, where `a`, `b`, and `c` are arrays, returns the values of `vcat(a, b, c)[ind]`.

This adapts to multi-dimensional arrays, respectively.

This is not exported.
"""
function get_lazy_vcat_indices(arr::AbstractVector, inds)
    nd = ndims(arr[1])
    if nd == 1
        outarr = Vector{eltype(arr[1])}(undef, length(inds))
    else
        outarr = Matrix{eltype(arr[1])}(undef, size(arr[1], 1), length(inds))
    end

    @inbounds for (i, ind) ∈ enumerate(inds)
        isfound = false
        for a ∈ arr
            n = nd == 1 ? length(a) : size(a, 2)
            if ind ≤ n
                if ndims(a) == 1
                    outarr[i] = a[ind]
                else
                    @views outarr[:,i] = a[:,ind]
                end
                isfound = true
                break
            end

            ind -= n
        end

        if !isfound
            throw(BoundsError)
        end
    end

    return outarr
end
