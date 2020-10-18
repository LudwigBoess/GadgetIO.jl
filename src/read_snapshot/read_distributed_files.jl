using ProgressMeter

function filter_positions(snap_file::String, corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real}, 
                          parttype::Integer)

    pos = read_snap(snap_file, "POS", parttype)
    sel = findall( ( corner_lowerleft[1] .<= pos[:,1] .< corner_upperright[1] ) .&
                   ( corner_lowerleft[2] .<= pos[:,2] .< corner_upperright[2] ) .&
                   ( corner_lowerleft[3] .<= pos[:,3] .< corner_upperright[3] ) )
    return sel
end

"""
    get_npart_to_read( snap_base::String, filter_function::Function)

Finds the number of particles in a snapshot directory that pass the `filter_function`.
"""
function get_npart_to_read( snap_base::String, filter_function::Function)

    # read the header of the zero'th file 
    h = read_header(select_file(snap_base, 0))

    # number of particles that fulfill the filter criteria
    N_part = 0

    for sub_snap = 0:(h.num_files-1)
        snap_file = select_file(snap_base, sub_snap)
        sel       = filter_function(snap_file)
        N_this_file = length(sel)
        @info "sub-snap $sub_snap: $N_this_file particles"
        N_part   += N_this_file
    end

    @info "Need to read $N_part particles"

    return N_part

end

"""
    read_blocks_over_all_files( snap_base::String, blocks::Array{String}, filter_function::Function; 
                                N_to_read::Integer=-1, parttype::Integer=0 )

Reads the specified blocks from all distributed files where particles pass the `filter_function`.
Per default the functions checks the number of particles that pass the filter. 
If that number is known in advance it can be given via the `N_to_read` keyword argument.
"""
function read_blocks_over_all_files(snap_base::String, blocks::Array{String}, filter_function::Function; 
                                    N_to_read::Integer=-1, parttype::Integer=0)

    h           = read_header(select_file(snap_base, 0))
    snap_info   = read_info(select_file(snap_base, 0))

    # default behaviour is to find the particles to read.
    # If this number is known in advance it can be given as an input parameter.
    if N_to_read == -1
        # check how many particles fulfill the filter criterion
        # -> needed for array pre-allocation to speed up read-in.
        N_to_read = get_npart_to_read(snap_base, filter_function)
    end

    # pre-allocate all data arrays in a dictionary
    d = allocate_data_dict(blocks, N_to_read, snap_info, false)

    # store the number of particles that have been read
    N_read = 0

    @showprogress "Reading..." for sub_snap = 0:(h.num_files-1)

        # select current file
        snap_file = select_file(snap_base, sub_snap)

        # find the particles in the file that fulfil the criterion
        sel         = filter_function(snap_file)
        N_this_file = length(sel)

        # read the blocks
        for block ∈ blocks

            # store the block in a dummy array
            dummy = read_snap(snap_file, block, parttype)
            
            # save the relevant particles in the dictionary
            for i = 1:N_this_file
                d[block][N_read+i,:] = dummy[sel[i][1],:]
            end

        end

        # count the number of particles already read
        N_read += N_this_file
    
    end

    return d
end