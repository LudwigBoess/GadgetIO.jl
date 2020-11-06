using ProgressMeter
using Base.Threads

function filter_positions(snap_file::String, corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real}, 
                          parttype::Integer)

    pos = read_snap(snap_file, "POS", parttype)
    sel = findall( ( corner_lowerleft[1] .<= pos[:,1] .< corner_upperright[1] ) .&
                   ( corner_lowerleft[2] .<= pos[:,2] .< corner_upperright[2] ) .&
                   ( corner_lowerleft[3] .<= pos[:,3] .< corner_upperright[3] ) )
    return sel
end



"""
    read_blocks_over_all_files( snap_base::String, blocks::Array{String}, filter_function::Function; 
                                N_to_read::Integer=-1, parttype::Integer=0 )

Reads the specified blocks from all distributed files where particles pass the `filter_function`.
Per default the functions checks the number of particles that pass the filter. 
If that number is known in advance it can be given via the `N_to_read` keyword argument.
"""
function read_blocks_over_all_files(snap_base::String, blocks::Array{String}, filter_function::Function; 
                                    read_positions=nothing, parttype::Integer=0, verbose::Bool=true )

    h           = read_header(select_file(snap_base, 0))
    snap_info   = read_info(select_file(snap_base, 0))

    # default behaviour is to find the particles to read.
    # If this number is known in advance it can be given as an input parameter.
    if isnothing(read_positions)
        # check how many particles fulfill the filter criterion and find the positions in the blocks.
        # -> needed for array pre-allocation to speed up read-in.
        read_positions = find_read_positions(snap_base, filter_function, verbose=verbose)
    end

    if verbose
        @info "Allocating storage dictionary..."
        t1 = time_ns()
    end
    # pre-allocate all data arrays in a dictionary
    d = allocate_data_dict(blocks, read_positions["N_part"], snap_info, false)
    
    if verbose
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    # store the number of particles that have been read
    N_read = 0

    key_entries = collect(keys(read_positions))
    files = key_entries[ key_entries .!= "N_part"]

    if verbose
        @info "Reading $(h.num_files) snapshots..."
        t1 = time_ns()
    end

    for file ∈ files

        # select current file
        filename = select_file(snap_base, file ) #parse(Int,file)

        # read header block
        h = read_header(filename)

        # read info block
        snap_info = read_info(filename)

        # find the location of the blocks in the current file
        file_block_positions = get_block_positions(filename)

        # read the blocks
        @threads for block ∈ blocks

            # read the info ot the curent block
            block_info = snap_info[getfield.(snap_info, :block_name) .== block][1]

            # add offset of particle types that should not be read
            offset = 0
            for i=1:size(h.npart)[1]
                if (block_info.is_present[i] > 0) & (h.npart[i] > 0) & ( i < parttype + 1)
                    offset += h.npart[i]
                end
            end

            read_block_with_offset!(d[block], N_read, filename, 
                        file_block_positions[block],
                        block_info, offset, 
                        read_positions[file]["index"],
                        read_positions[file]["n_to_read"])

        end
        # count the number of particles already read
        N_read += sum(read_positions[file]["n_to_read"])

        if verbose
            @info "Read $N_read / $(read_positions["N_part"]) particles"
        end
        
    
    end

    if verbose
        t2 = time_ns()
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    return d
end