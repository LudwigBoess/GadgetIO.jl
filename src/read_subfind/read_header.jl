"""
    struct SubfindHeader

Contains the data of the `HEAD` block in the subfind output

# Fields
| Name                                 | Meaning                                                                                |
| :----------------------------------  | :------------------------------------------------------------------------------------- |
| `nhalos::Int32`                      | number of halos in the output file                                                     |
| `nsubhalos::Int32`                   | number of subhalos in the output file                                                  |
| `nfof::Int32`                        | number of particles in the FoF                                                         |
| `ngroups::Int32`                     | number of large groups in the output file                                              |
| `time::Float64`                      | time / scale factor of the simulation                                                  |
| `z::Float64`                         | redshift of the simulation                                                             |
| `tothalos::UInt32`                   | total number of halos over all output files                                            |
| `totsubhalos::UInt32`                | total number of subhalos over all output files                                         |
| `totfof::UInt32`                     | total number of particles in the FoF                                                   |
| `totgroups::UInt32`                  | 1 if simulation was run with cooling, else 0                                           |
| `num_files::Int32`                   | number of files over which subfind data is distributed                                 |
| `boxsize::Float64`                   | total size of the simulation box                                                       |
| `omega_0::Float64`                   | Omega matter                                                                           |
| `omega_l::Float64`                   | Omega dark enery                                                                       |
| `h0::Float64`                        | little h                                                                               |
| `flag_doubleprecision::Int32`        | 1 if snapshot is in double precision, else 0                                           |
| `flag_ic_info::Int32`                | 1 if initial snapshot file contains an info block, else 0                              |

"""
struct SubfindHeader <: AbstractGadgetHeader
    nhalos::Int32                       # number of halos in the output file
    nsubhalos::Int32                    # number of subhalos in the output file
    nfof::Int32                         # number of particles in the FoF
    ngroups::Int32                      # number of large groups in the output file
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    tothalos::UInt32                    # total number of halos over all output files
    totsubhalos::UInt32                 # total number of subhalos over all output files
    totfof::UInt32                      # total number of particles in the FoF
    totgroups::UInt32                   # total number of large groups over all output files
    num_files::Int32                    # number of files
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32
end


"""
    read_subfind_header(filename::String)

Reads the header of a subfind file or file base (without .0, .1, etc.) into a [`SubfindHeader`](@ref) struct.
"""
function read_subfind_header(filename::String)

    # read the SnapshotHeader
    h = read_header(filename)

    # convert it to a SubfindHeader
    convert_header(h)
end


"""
    convert_header(h::SnapshotHeader)

Converts a [`SnapshotHeader`](@ref) to a [`SubfindHeader`](@ref).
"""
function convert_header(h::SnapshotHeader)

    # construct total number of halos
    tothalos    = get_total_particles(h.nall[1],    h.npartTotalHighWord[1])
    totsubhalos = get_total_particles(h.nall[2],    h.npartTotalHighWord[2])
    totfof      = get_total_particles(h.nall[3],    h.npartTotalHighWord[3])
    totgroups   = get_total_particles(h.nall[4],    h.npartTotalHighWord[4])

    return SubfindHeader(h.npart[1], h.npart[2], h.npart[3], h.npart[4],
                         h.time, h.z,
                         tothalos, totsubhalos, totfof, totgroups,
                         h.num_files,
                         h.boxsize, h.omega_0, h.omega_l, h.h0,
                         h.flag_doubleprecision,
                         h.flag_ic_info)

end