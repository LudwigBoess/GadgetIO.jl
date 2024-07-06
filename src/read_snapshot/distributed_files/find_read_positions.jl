using ProgressMeter

"""
    reduce_read_positions(sel::Array{<:Integer})

Reduces the individual read positions by finding sub-ranges. 
Returns an array of read indices and an array with the number of entries per index.
"""
function reduce_read_positions(sel::Array{<:Integer})

    Npart     = size(sel,1)
    index     = zeros(Int64, Npart) 
    n_to_read = zeros(Int64, Npart)

    index[1]     = sel[1] - 1
    n_to_read[1] = 1
    Nentries     = 1

    @inbounds for  i=2:Npart
        if (sel[i] == index[Nentries] + n_to_read[Nentries] + 1)
            n_to_read[Nentries] += 1
        else
            Nentries            += 1
            index[Nentries]     = sel[i] - 1
            n_to_read[Nentries] += 1
        end
    end

    resize!(index, Nentries)
    resize!(n_to_read, Nentries)

    return index, n_to_read
end


"""
    find_read_positions( snap_base::String, filter_function::Function)

Finds the number of particles and their storage location in a snapshot directory that pass the `filter_function`.
"""
function find_read_positions( snap_base::String, filter_function::Function;
                              verbose::Bool=true)

    if verbose
        @info "Getting number of particles to read..."
        t1 = time_ns()
    end

    # get number of sub-snapshots (if any) 
    snap_file = select_file(snap_base, 0)
    if HDF5.ishdf5(snap_file)
        num_files = h5readattr(snap_file, "Header")["NumFilesPerSnapshot"]
    else
        num_files = read_header(snap_file).num_files
    end

    # number of particles that fulfill the filter criteria
    N_part = 0

    # storage dictionary
    d = Dict()

    for sub_snap = 0:(num_files-1)
        snap_file   = select_file(snap_base, sub_snap)
        sel         = filter_function(snap_file)
        N_this_file = size(sel,1)

        # read number of particles in file
        if HDF5.ishdf5(snap_file)
            npart = h5readattr(snap_file, "Header")["NumPart_ThisFile"]
        else
            npart = read_header(snap_file).npart
        end

        # store read positions of particles are in the file
        if N_this_file > 0
            index, n_to_read = reduce_read_positions(sel)

            d[sub_snap] = Dict( "index" => index, "n_to_read" => n_to_read)
        end

        if verbose
            @info "sub-snap $sub_snap: $N_this_file particles"
        end
        N_part   += N_this_file
    end

    if verbose
        t2 = time_ns()
        snap_file   = select_file(snap_base, 0)
        @info "Need to read $N_part particles"
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    d["N_part"] = N_part

    return d

end


"""
    find_read_positions( snap_base::String, geometry::AbstractGadgetGeometry;
                         parttype::Integer=0,
                         verbose::Bool=true )

Finds the number of particles and their storage location in a snapshot directory that are contained in the space defined by `geometry`.
"""
function find_read_positions( snap_base::String, geometry::AbstractGadgetGeometry;
                              parttype::Integer=0,
                              verbose::Bool=true)

    if verbose
        @info "Getting number of particles to read..."
        t1 = time_ns()
    end

    # get number of sub-snapshots (if any) 
    snap_file = select_file(snap_base, 0)
    if HDF5.ishdf5(snap_file)
        num_files = h5readattr(snap_file, "Header")["NumFilesPerSnapshot"]
    else
        num_files = read_header(snap_file).num_files
    end

    # number of particles that fulfill the filter criteria
    N_part = 0

    # storage dictionary
    d = Dict()

    for sub_snap = 0:(num_files-1)
        snap_file = select_file(snap_base, sub_snap)

        # read the position block
        if HDF5.ishdf5(snap_file)
            pos   = h5read(snap_file, "/PartType$parttype")["Coordinates"]
        else
            pos   = read_block(snap_file, "POS"; parttype)
        end

        # find positions of particles that are contained in the geometry
        sel         = get_geometry_mask(geometry, pos)
        N_this_file = size(sel,1)

        # read number of particles in file
        if HDF5.ishdf5(snap_file)
            npart = h5readattr(snap_file, "Header")["NumPart_ThisFile"]
        else
            npart = read_header(snap_file).npart
        end

        # store read positions of particles are in the file
        if N_this_file > 0
            index, n_to_read = reduce_read_positions(sel)

            d[sub_snap] = Dict( "index" => index, "n_to_read" => n_to_read)
        end

        if verbose
            @info "sub-snap $sub_snap: $N_this_file / $(npart[parttype+1]) particles"
        end
        N_part   += N_this_file
    end

    if verbose
        t2 = time_ns()
        snap_file   = select_file(snap_base, 0)
        if HDF5.ishdf5(snap_file)
            N_all = Int64(h5readattr(snap_file, "Header")["NumPart_Total"][parttype+1])
        else
            h = read_header(snap_file)
            N_all = Int64(get_total_particles(h, parttype))
        end
        N_per = N_part / N_all * 100.0
        @info "Need to read $N_part / $N_all particles -> $(@sprintf("%0.2f", N_per)) %"
        @info "  elapsed: $(output_time(t1,t2)) s"
    end

    d["N_part"] = N_part

    return d

end

"""
    save_read_positions(save_file::String, read_positions::Dict)

Saves the relevant read-in positions to a binary file.
"""
function save_read_positions(save_file::String, read_positions::Dict)

    f = open(save_file, "w")
    N_part = read_positions["N_part"]

    delete!(read_positions, "N_part")
    files = keys(read_positions)
    N_files = length(keys(read_positions))

    write(f, Int64(N_part))
    write(f, Int64(N_files))

    for file ∈ files
        N_entries = size(read_positions[file]["index"],1)
        write(f, Int64(file))
        write(f, Int64(N_entries))
        write(f, Int64.(read_positions[file]["index"]))
        write(f, Int64.(read_positions[file]["n_to_read"]))
    end

    close(f)

end

"""
    load_read_positions(save_file::String)

Loads the relevant read-in positions from a binary file.
"""
function load_read_positions(save_file::String)
    
    read_positions = Dict()

    f = open(save_file, "r")

    read_positions["N_part"] = read(f, Int64)

    N_files = read(f, Int64)

    for _ ∈ 1:N_files

        file      = read(f, Int64)
        N_entries = read(f, Int64)
        index     = read!(f, Array{Int64,1}(undef, N_entries))
        n_to_read = read!(f, Array{Int64,1}(undef, N_entries))

        read_positions[file] = Dict("index" => index, "n_to_read" => n_to_read)

    end

    close(f)

    return read_positions
end
