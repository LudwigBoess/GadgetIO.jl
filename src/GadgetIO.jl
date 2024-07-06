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
include(joinpath("read_snapshot", "utility", "default_info_lines.jl"))
include(joinpath("read_snapshot", "utility", "allocate.jl"))
include(joinpath("read_snapshot", "utility", "snapshot_utilities.jl"))
include(joinpath("read_snapshot", "utility", "block_position.jl"))
include(joinpath("read_snapshot", "utility", "total_particles.jl"))
include(joinpath("read_snapshot", "read_block.jl"))
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


"""
    Precompile steps
"""

using PrecompileTools    # this is a small dependency
using Downloads

@setup_workload begin

    data_path = joinpath(@__DIR__, "..", "precompile_data")
    sub_base = joinpath(data_path, "sub_002")
    snap_base = joinpath(data_path, "snap_002")
    key_file = joinpath(data_path, "snap_144.0.key")
    key_index = joinpath(data_path, "snap_144.key.index")

    low_bounds, high_bounds = [0, 4, 7, 10], [1, 5, 8, 16]
    x0 = [-100.0, -100.0, -100.0]
    x1 = [1_000.0, 1_000.0, 1_000.0]

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)

        # header and info 
        h = read_header(snap_base)
        info = read_info(snap_base * ".0")

        # blocks 
        pos = read_block(snap_base * ".0", "POS", parttype=0)
        mass = read_block(snap_base * ".0", "MASS", parttype=0)
        pos = read_block(snap_base * ".0", "POS", parttype=1)
        mass = read_block(snap_base * ".0", "MASS", parttype=1)

        # particles in box
        center = Float32[0.5, 0.5, 0.5]
        rvir = 0.3f0
        pos = read_particles_in_volume(snap_base, "POS", center, rvir, use_keys=false, parttype=1, verbose=false)

        # particles in geometry
        cube = GadgetCube(center .- rvir, center .+ rvir)
        pos = read_particles_in_geometry(snap_base, "POS", cube, use_keys=false, parttype=1, verbose=false)
        sphere = GadgetSphere(center, rvir)
        data = read_particles_in_geometry(snap_base, "POS", sphere, use_keys=false, parttype=1, verbose=false)

        # read positions 
        ff(filename) = filter_cube(filename, center .- rvir, center .+ rvir, parttype=1)
        read_positions = find_read_positions(snap_base, ff; verbose=false)
        data = read_blocks_filtered(snap_base, ["POS", "MASS"], verbose=false; read_positions)

        # filter function 
        data = read_blocks_filtered(snap_base, ["POS", "MASS"], filter_function=ff, verbose=false)

        # subfind 
        mtop = read_subfind(sub_base, "MTOP")
        mtop2, haloids = read_subfind(sub_base, "MTOP", return_haloid=true)

        # PH-keys 
        low_list, high_list, file_list = GadgetIO.read_key_index(key_index)
        h_key = GadgetIO.read_keyheader(key_file)

        GadgetIO.peano_hilbert_key(h_key.bits, 0, 0, 0)
        GadgetIO.peano_hilbert_key(h_key.bits, 0, 0, 1)
        GadgetIO.peano_hilbert_key(h_key.bits, 0, 1, 1)

        GadgetIO.get_index_bounds([4, 5, 8, 12, 15], low_bounds, high_bounds)

        keylist = GadgetIO.get_keylist(h_key, x0, x1)

        GadgetIO.get_int_pos(1000.5, h_key.domain_corners[1], h_key.domain_fac, h_key.bits)
    end

end

end # module
