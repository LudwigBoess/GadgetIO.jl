__precompile__()

module GadgetIO

    # functions to read snapshots
    include(joinpath(dirname(@__FILE__), "read_snapshot", "gadget_types.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_header.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "snapshot_utilities.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_format_1.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_format_2.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_snapshot.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_subfind.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_particle_in_box.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_particles_in_halo.jl"))

    # functions to write snapshots
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))

    export Header, Info_Line,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           read_info,
           block_present,
           read_snap,
           read_block_by_name,      # similar to readnew.pro by Klaus Dolag
           read_header,
           read_particles_in_box,
           read_particles_in_volume,
           read_particles_in_halo,
           get_index_list,

           # subfind read
           HaloID,
           read_subfind_header,
           read_subfind,
           find_most_massive_halo,
           filter_subfind,
           read_ids_in_halo,
           read_particles_by_id,

           # write snapshot functions
           write_header,
           write_block

end # module
