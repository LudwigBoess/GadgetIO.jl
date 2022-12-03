"""
    read_halo_prop(filebase, haloid::HaloID, blockname::AbstractString; verbose::Bool=true)

Get halo property from block `blockname` by halo id. Returns an array or scalar depending on the block type.
"""
function read_halo_prop(filebase, blockname::AbstractString, haloid::HaloID; verbose::Bool=true)
    
    if verbose
        @info "Reading property $blockname of halo $haloid"
    end

    # get full specified block for all halos in file
    sub_input = select_file(filebase, haloid.file)

    prop = read_subfind(sub_input, blockname, offset=haloid.id-1, n_to_read=1)

    # prop is either a length=1 vector or a dim×1 matrix, so it is extracted to scalar or vector
    if ndims(prop) == 1
        return prop[1]
    else
        return vec(prop)
    end
end

"""
    read_halo_prop(filebase, i_global::Integer, blockname::AbstractString; verbose::Bool=true)

Get halo property from block `blockname` by global halo index `i_global` (zero-based index).
"""
function read_halo_prop(filebase, blockname::AbstractString, i_global::Integer; verbose::Bool=true)
    val, id = read_halo_prop_and_id(filebase, i_global, blockname; verbose)
    return val
end


"""
    read_halo_prop_and_id(filebase, i_global::Integer, blockname::AbstractString; kwargs...)

Get halo property and [`HaloID`](@ref) from block `blockname` by global halo index `i_global` (zero-based index).
`nfiles` should generally be set to `h.num_files`, obtained from `read_header`.
When nfiles is not passed, it is read automatically from the header.
"""
function read_halo_prop_and_id(sub_base, blockname::AbstractString, i_global::Integer; 
                                verbose::Bool=true)
    if verbose
        @info "Reading property $blockname of halo at index $i_global"
    end

    # blocks are type specific so we can use this to make our life easier
    parttype = subfind_block_parttype(sub_base, blockname)

    # convert global halo idx to HaloID
    haloid = global_idxs_to_halo_id(sub_base, [i_global]; parttype)[1]


    return read_halo_prop(sub_base, blockname, haloid; verbose), haloid
end
