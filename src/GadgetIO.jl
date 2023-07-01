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


"""
    Precompile steps
"""

using PrecompileTools    # this is a small dependency
using Downloads

@setup_workload begin

    # download the snapshot data for precompilation
    @info "downloading data..."

    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.0", "./sub_002.0")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.1", "./sub_002.1")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.2", "./sub_002.2")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/sub_002.3", "./sub_002.3")

    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.0", "./snap_002.0")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.1", "./snap_002.1")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.2", "./snap_002.2")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.3", "./snap_002.3")

    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.0.key", "./snap_002.0.key")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.1.key", "./snap_002.1.key")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.2.key", "./snap_002.2.key")
    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.3.key", "./snap_002.3.key")

    Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_002.key.index", "./snap_002.key.index")


    @info "done!"

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)

        # header and info 
        h = read_header("snap_002")
        info = read_info("snap_002.0")

        # blocks 
        pos = read_block("snap_002.0", "POS", parttype=0)
        mass = read_block("snap_002.0", "MASS", parttype=0)
        mass_full = read_block("snap_002.0", "MASS", parttype=-1)


        # particles in box
        center = Float32[3978.9688, -95.40625, -8845.25]
        rvir = 118.76352
        pos = read_particles_in_volume("snap_002", "POS", center, rvir, use_keys=false, parttype=1)
        pos = read_particles_in_volume("snap_002", "POS", center, rvir, use_keys=true, parttype=1, verbose=true)
        id = read_particles_in_volume("snap_002", "ID", center, rvir, use_keys=true, parttype=1, verbose=false)

        # particles in geometry
        center = [3978.9688, -95.40625, -8845.25]
        rvir = 118.76352
        cube = GadgetCube(center .- rvir, center .+ rvir)
        pos = read_particles_in_geometry("snap_002", "POS", cube, use_keys=false, parttype=1)
        sphere = GadgetSphere(center, rvir)
        data = read_particles_in_geometry("snap_002", "POS", sphere, use_keys=false, parttype=1)

        # particles in halo
        pos = read_particles_in_halo("snap_002", "POS", "sub_002", HaloID(0, 4), use_keys=false)
        ids = UInt32[0x000028fc, 0x00002594, 0x00002963, 0x00002681, 0x00001af4, 0x00001ff1, 0x000022d7, 0x00002267, 0x000029c0, 0x0000277b]
        pos = read_particles_by_id("snap_002", ids, "POS")

        # read positions 
        ff(filename) = filter_cube(filename, center .- rvir, center .+ rvir, parttype=1)
        read_positions = find_read_positions("snap_002", ff)

        # subfind 
        center, rvir, haloid = find_most_massive_halo("sub_002", 4)
        mtop = read_subfind("sub_002", "MTOP")
        mtop2, haloids = read_subfind("sub_002", "MTOP", return_haloid=true)
    end

    rm("sub_002.0")
    rm("sub_002.1")
    rm("sub_002.2")
    rm("sub_002.3")
    rm("snap_002.0")
    rm("snap_002.1")
    rm("snap_002.2")
    rm("snap_002.3")
    rm("snap_002.0.key")
    rm("snap_002.1.key")
    rm("snap_002.2.key")
    rm("snap_002.3.key")
    rm("snap_002.key.index")

end

end # module
