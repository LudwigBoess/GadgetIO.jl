module GadgetIO

using Printf

function output_time(t1, t2)
    return @sprintf("%0.3e", Float64((t2 - t1)) * 1.e-9)
end

"""
    AbstractGadgetHeader

Abstract type to unify snapshot and subfind headers.
"""
abstract type AbstractGadgetHeader end

# functions to read snapshots
include(joinpath("read_snapshot", "utility", "domain_slices.jl"))
include(joinpath("read_snapshot", "read_header.jl"))
include(joinpath("read_snapshot", "read_info.jl"))
include(joinpath("read_snapshot", "utility", "allocate.jl"))
include(joinpath("read_snapshot", "utility", "snapshot_utilities.jl"))
include(joinpath("read_snapshot", "utility", "block_position.jl"))
include(joinpath("read_snapshot", "utility", "total_particles.jl"))
include(joinpath("read_snapshot", "utility", "mass_assignment.jl"))
include(joinpath("read_snapshot", "read_format_1.jl"))
include(joinpath("read_snapshot", "read_format_2.jl"))
include(joinpath("read_snapshot", "read_snapshot.jl"))

include(joinpath("read_subfind", "read_header.jl"))
include(joinpath("read_subfind", "halo_ids.jl"))
include(joinpath("read_subfind", "read_subfind.jl"))
include(joinpath("read_subfind", "read_halo_prop.jl"))
include(joinpath("read_subfind", "filter_subfind.jl"))

include(joinpath("read_snapshot", "particles_in_box", "read_key_files.jl"))
include(joinpath("read_snapshot", "particles_in_box", "peano_hilbert.jl"))
include(joinpath("read_snapshot", "particles_in_box", "utility.jl"))
include(joinpath("read_snapshot", "particles_in_box", "read_particle_in_box.jl"))
include(joinpath("read_snapshot", "particles_in_geometry", "abstract_geometry.jl"))
include(joinpath("read_snapshot", "particles_in_geometry", "cylinder.jl"))
include(joinpath("read_snapshot", "particles_in_geometry", "sphere.jl"))
include(joinpath("read_snapshot", "particles_in_geometry", "shell.jl"))
include(joinpath("read_snapshot", "particles_in_geometry", "cube.jl"))
include(joinpath("read_snapshot", "particles_in_geometry", "read_particles_in_geometry.jl"))
include(joinpath("read_snapshot", "particles_in_halo", "read_particles_by_id.jl"))
include(joinpath("read_snapshot", "particles_in_halo", "read_particles_in_halo.jl"))
include(joinpath("read_snapshot", "distributed_files", "find_read_positions.jl"))
include(joinpath("read_snapshot", "distributed_files", "read_distributed_files.jl"))

include(joinpath("timer_outputs", "read_balance.jl"))

# functions to write snapshots
include(joinpath("write_snapshot", "write_snap.jl"))

export AbstractGadgetHeader,
    SnapshotHeader, SubfindHeader,
    InfoLine,       # types
    head_to_dict,
    snap_to_dict,
    head_to_struct,
    print_blocks,
    read_info,
    block_present,
    read_snap,
    read_block,      # similar to readnew.pro by Klaus Dolag
    read_header,

    # large simulations
    read_particles_in_box,
    read_particles_in_volume,
    read_particles_in_geometry,
    read_particles_in_halo,
    get_index_list,
    get_npart_to_read,
    filter_ids,
    read_blocks_filtered,
    find_read_positions,
    save_read_positions,
    load_read_positions,
    get_total_particles,

    # geometries
    AbstractGadgetGeometry,
    GadgetCylinder,
    GadgetSphere,
    GadgetShell,
    GadgetCube,

    # domain slices
    filter_cube,
    filter_cylinder,

    # subfind read
    HaloID,
    SubfindHeader,
    read_subfind_header,
    read_subfind,
    find_most_massive_halo,
    filter_subfind,
    save_halo_ids,
    load_halo_ids,
    read_ids_in_halo,
    read_particles_by_id,
    read_halo_prop,
    read_halo_prop_and_id,
    global_idxs_to_halo_id,
    halo_ids_to_read_positions,

    # write snapshot functions
    write_header,
    write_block,
    write_info_block,

    # timer files 
    parse_balance,
    print_performance

end # module
