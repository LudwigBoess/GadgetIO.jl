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
        dim   = snap_info[getfield.(snap_info, :block_name) .== block][1].n_dim
        dtype = snap_info[getfield.(snap_info, :block_name) .== block][1].data_type
        if dim == 1
            d[block] = Vector{dtype}(undef, N_to_read)
        else
            d[block] = Matrix{dtype}(undef, dim, N_to_read)
        end
    end
    
    # allocate mass array, if it's not in a block
    if no_mass_block
        d["MASS"] = Vector{Float32}(undef, N_to_read)
    end

    return d
end
