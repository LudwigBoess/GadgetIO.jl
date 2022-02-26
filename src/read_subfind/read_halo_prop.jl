"""
    get_prop_from_block(block::Array, index::Integer)

Get the full property from `block` at `index`. Returns an array or scalar depending on the block type.
"""
function get_prop_from_block(block::Array, index::Integer)
    return ndims(block) == 1 ? block[index] : block[:, index]
end



"""
    read_halo_prop(filebase, haloid::HaloID, blockname::AbstractString; verbose::Bool=true)

Get halo property from block `blockname` by halo id. Returns an array or scalar depending on the block type.
"""
function read_halo_prop(filebase, haloid::HaloID, blockname::AbstractString; verbose::Bool=true)
    
    if verbose
        @info "Reading property $blockname of halo $haloid"
    end

    # get full specified block for all halos in file
    sub_input = select_file(filebase, haloid.file)

    prop = read_subfind(sub_input, blockname; offset=haloid.id-1, nread=1)

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
function read_halo_prop(filebase, i_global::Integer, blockname::AbstractString; verbose::Bool=true)
    val, id = read_halo_prop_and_id(filebase, i_global, blockname; verbose)
    return val
end

"""
    read_halo_prop_and_id(filebase, i_global::Integer, blockname::AbstractString; kwargs...)
    read_halo_prop_and_id(filebase, i_global::Integer, blockname::AbstractString, nfiles::Integer; verbose::Bool=true)

Get halo property and [`HaloID`](@ref) from block `blockname` by global halo index `i_global` (zero-based index).
`nfiles` should generally be set to `h.num_files`, obtained from `read_header`.
When nfiles is not passed, it is read automatically from the header.
"""
function read_halo_prop_and_id(filebase, i_global::Integer, blockname::AbstractString; kwargs...)
    h = read_subfind_header(filebase)
    nfiles = h.num_files

    read_halo_prop_and_id(filebase, i_global, blockname, nfiles; kwargs...)
end

function read_halo_prop_and_id(filebase, i_global::Integer, blockname::AbstractString, nfiles::Integer; verbose::Bool=true)
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

        len = read_subfind_length(sub_input, blockname)
        if len ≤ i_global # halo is not in current file
            i_global -= len
        else
            prop = read_subfind(sub_input, blockname; offset=i_global, nread=1)
            haloid = HaloID(i, i_global + 1)

            # prop is either a length=1 vector or a dim×1 matrix, so it is extracted to scalar or vector
            if ndims(prop) == 1
                return prop[1], haloid
            else
                return vec(prop), haloid
            end
        end
    end
    error("Halo at index $i_global0 does not exist!")
end
