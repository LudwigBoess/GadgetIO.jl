"""
    get_requested_info(snap_info::Array{InfoLine}, block::AbstractString)

Checks if the `InfoLine` is present for the requested block, or if the `MASS` block is requested, but not in the `INFO` block. 
"""
function get_requested_info(snap_info::Array{InfoLine}, block::AbstractString)

     # read the info ot the curent block
     if !isempty(snap_info[getfield.(snap_info, :block_name) .== block])
        return snap_info[getfield.(snap_info, :block_name) .== block][1]
    else
        if block == "MASS"
            # ! TODO: more stable form of this!
            return InfoLine("MASS", Float32, 1, [1, 1, 1, 1, 1, 1])
        else
            error("INFO missing for requested block: $block !")
        end
    end

end

"""
    allocate_data_dict( blocks::Array{String}, N_to_read::Integer, 
                        snap_info::Array{InfoLine}, no_mass_block::Bool )

Helper function to allocate the data `Dict`.
"""
function allocate_data_dict(blocks::Array{String}, N_to_read::Integer, 
                            snap_info::Array{InfoLine}, no_mass_block::Bool)

    # prepare dictionary for particle storage
    d = Dict{String, VecOrMat{T} where T}()

    for block ∈ blocks
        
        block_info = get_requested_info(snap_info, block)

        if block_info.n_dim == 1
            d[block] = Vector{block_info.data_type}(undef, N_to_read)
        else
            d[block] = Matrix{block_info.data_type}(undef, block_info.n_dim, N_to_read)
        end
    end

    if no_mass_block
        d["MASS"] = Vector{Float32}(undef, N_to_read)
    end

    return d
end

"""
    allocate_data_dict( blocks::Array{String}, N_to_read::Integer, 
                        snap_info::Array{InfoLine}, no_mass_block::Bool )

Helper function to allocate the data `Dict`.
"""
function allocate_data_dict(blocks::Array{String}, N_to_read::Integer, 
                            group::HDF5.Group, no_mass_block::Bool)

    # prepare dictionary for particle storage
    d = Dict{String, VecOrMat{T} where T}()

    for block ∈ blocks
        
        # if 1D allocate Vector
        if ndims(group[block]) == 1
            d[block] = Vector{eltype(group[block])}(undef, N_to_read)
        else
            ndim = size(group[block])[1]
            d[block] = Matrix{eltype(group[block])}(undef, ndim, N_to_read)
        end
    end

    if no_mass_block
        d["Masses"] = Vector{Float32}(undef, N_to_read)
    end

    return d
end

function allocate_data_array(info::InfoLine, N_to_read::Integer)

    if info.n_dim == 1
        return Vector{info.data_type}(undef, N_to_read)
    else
        return Matrix{info.data_type}(undef, info.n_dim, N_to_read)
    end
end



