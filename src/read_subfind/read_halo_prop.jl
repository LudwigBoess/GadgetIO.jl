"""
    get_prop_from_block(block::Array, index::Integer)
Get the full property from `block` at `index`. Returns an array or scalar depending on the block type.
"""
function get_prop_from_block(block::Array, index::Integer)
    return ndims(block) == 1 ? block[index] : block[:, index]
end



"""
    read_halo_prop(filebase::String, haloid::HaloID, blockname::String; verbose::Bool=true)
Get halo property from block `blockname` by halo id. Returns an array or scalar depending on the block type.
"""
function read_halo_prop(filebase::String, haloid::HaloID, blockname::String; verbose::Bool=true)
    
    if verbose
        @info "Reading property $blockname of halo $haloid"
    end

    # get full specified block for all halos in file
    sub_input = select_file(filebase, haloid.file)
    block = read_subfind(sub_input, blockname)

    # return array or scalar for given halo
    return get_prop_from_block(block, haloid.id)
end


"""
    read_halo_prop_and_id(filebase::String, i_global::Integer, blockname::String, nfiles::Integer=1; verbose::Bool=true)
Get halo property from block `blockname` by global halo index `i_global`.
`nfiles` should generally be set to `h.num_files`, obtained from `read_header`.
"""
function read_halo_prop_and_id(filebase::String, i_global::Integer, blockname::String, nfiles::Integer=1; verbose::Bool=true)
    if verbose
        @info "Reading property $blockname of halo at index $i_global"
    end

    # store here for error handling
    i_global0 = i_global

    # run through files until halo at index i_global is reached
    for i in 0:nfiles - 1
        sub_input = select_file(filebase, i)

        if verbose
            @info "Reading file $i of $(nfiles - 1)"
        end

        block = read_subfind(sub_input, blockname)
        len = size(block, 1)
        if len ≤ i_global # halo is not in current file
            i_global -= len
        else
            return get_prop_from_block(block, i_global + 1), HaloID(i, i_global + 1)
        end
    end
    error("Halo at index $i_global0 does not exist!")
end