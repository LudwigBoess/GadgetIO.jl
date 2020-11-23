__precompile__()

module GadgetIO

    using Printf

    function output_time(t1, t2)
        return @sprintf("%0.3e", Float64((t2-t1))*1.e-9)
    end
    
    # functions to read snapshots
    include(joinpath("read_snapshot", "read_header.jl"))
    include(joinpath("read_snapshot", "read_info.jl"))
    include(joinpath("read_snapshot", "utility", "allocate_dict.jl"))
    include(joinpath("read_snapshot", "utility", "snapshot_utilities.jl"))
    include(joinpath("read_snapshot", "read_format_1.jl"))
    include(joinpath("read_snapshot", "read_format_2.jl"))
    include(joinpath("read_snapshot", "read_snapshot.jl"))
    include(joinpath("read_subfind",  "read_subfind.jl"))
    include(joinpath("read_subfind",  "read_halo_prop.jl"))
    include(joinpath("read_subfind",  "filter_subfind.jl"))
    include(joinpath("read_snapshot", "distributed_files", "find_read_positions.jl"))
    include(joinpath("read_snapshot", "distributed_files", "read_distributed_files.jl"))
    include(joinpath("read_snapshot", "particles_in_box",  "read_key_files.jl"))
    include(joinpath("read_snapshot", "particles_in_box",  "peano_hilbert.jl"))
    include(joinpath("read_snapshot", "particles_in_box",  "utility.jl"))
    include(joinpath("read_snapshot", "particles_in_box",  "read_particle_in_box.jl"))
    include(joinpath("read_snapshot", "particles_in_halo", "read_particles_by_id.jl"))
    include(joinpath("read_snapshot", "particles_in_halo", "read_particles_in_halo.jl"))

    # functions to write snapshots
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))

    export SnapshotHeader, InfoLine,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           read_info,
           block_present,
           read_snap,
           read_block,      # similar to readnew.pro by Klaus Dolag
           read_header,

           # large simulations
           read_particles_in_box,
           read_particles_in_volume,
           read_particles_in_halo,
           get_index_list,
           get_npart_to_read,
           read_blocks_over_all_files,
           find_read_positions,
           save_read_positions,
           load_read_positions,
           filter_positions,


           # subfind read
           HaloID,
           read_subfind_header,
           read_subfind,
           find_most_massive_halo,
           filter_subfind,
           read_ids_in_halo,
           read_particles_by_id,
           read_halo_prop,
           read_halo_prop_and_id,

           # write snapshot functions
           write_header,
           write_block

end # module
