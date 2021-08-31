"""
    mutable struct SnapshotHeader( [ **Fields ])

Contains the data of the `HEAD` block of a Gadget snapshot.

# Fields
| Name                                 | Meaning                                                                                |
| :-------------------------------     | :------------------------------------------------------------------------------------- |
| `npart::Vector{Int32}`               | an array of particle numbers per type in this snapshot                                 |
| `massarr::Vector{Float64}`           | an array of particle masses per type in this snapshot - if zero: MASS block present    |
| `time::Float64`                      | time / scale factor of the simulation                                                  |
| `z::Float64`                         | redshift of the simulation                                                             |
| `flag_sfr::Int32`                    | 1 if simulation was run with star formation, else 0                                    |
| `flag_feedback::Int32`               | 1 if simulation was run with stellar feedback, else 0                                  |
| `nall::Vector{UInt32}`               | total number of particles in the simulation                                            |
| `flag_cooling::Int32`                | 1 if simulation was run with cooling, else 0                                           |
| `num_files::Int32`                   | number of snapshots over which the simulation is distributed                           |
| `omega_0::Float64`                   | Omega matter                                                                           |
| `boxsize::Float64`                   | total size of the simulation box                                                       |
| `omega_l::Float64`                   | Omega dark enery                                                                       |
| `h0::Float64`                        | little h                                                                               |
| `flag_stellarage::Int32`             | 1 if simulation was run with stellar age, else 0                                       |
| `flag_metals::Int32`                 | 1 if simulation was run with metals, else 0                                            |
| `npartTotalHighWord::Vector{UInt32}` | weird                                                                                  |
| `flag_entropy_instead_u::Int32`      | 1 if snapshot U field contains entropy instead of internal energy, else 0              |
| `flag_doubleprecision::Int32`        | 1 if snapshot is in double precision, else 0                                           |
| `flag_ic_info::Int32`                | 1 if initial snapshot file contains an info block, else 0                              |
| `lpt_scalingfactor::Float32`         | factor to use second order ic generation                                               |
| `fill::Vector{Int32}`                | the HEAD block needs to be filled with zeros to have a size of 256 bytes               |
"""
mutable struct SnapshotHeader <: AbstractGadgetHeader
    npart::Vector{Int32}                # an array of particle numbers per type in this snapshot
    massarr::Vector{Float64}            # an array of particle masses per type in this snapshot - if zero: MASS block present
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    flag_sfr::Int32                     # 1 if simulation was run with star formation, else 0
    flag_feedback::Int32                # 1 if simulation was run with stellar feedback, else 0
    nall::Vector{UInt32}                # total number of particles in the simulation
    flag_cooling::Int32                 # 1 if simulation was run with cooling, else 0
    num_files::Int32                    # number of snapshots over which the simulation is distributed
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_stellarage::Int32              # 1 if simulation was run with stellar age, else 0
    flag_metals::Int32                  # 1 if simulation was run with metals, else 0
    npartTotalHighWord::Vector{UInt32}  # weird
    flag_entropy_instead_u::Int32       # 1 if snapshot U field contains entropy instead of internal energy, else 0
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32                 # 1 if initial snapshot file contains an info block, else 0
    lpt_scalingfactor::Float32          # factor to use second order ic generation
    fill::Vector{Int32}                 # the HEAD block needs to be filled with zeros to have a size of 256 bytes

    function SnapshotHeader(npart::Vector{Int32}=Int32.([0,0,0,0,0,0]),
           massarr::Vector{Float64}=zeros(6),
           time::Float64=0.,
           z::Float64=0.,
           flag_sfr::Int32=Int32(0),
           flag_feedback::Int32=Int32(0),
           nall::Vector{UInt32}=UInt32.([0,0,0,0,0,0]),
           flag_cooling::Int32=Int32(0),
           num_files::Int32=Int32(0),
           boxsize::Float64=0.,
           omega_0::Float64=0.,
           omega_l::Float64=0.,
           h0::Float64=0.,
           flag_stellarage::Int32=Int32(0),
           flag_metals::Int32=Int32(0),
           npartTotalHighWord::Vector{UInt32}=UInt32.([0,0,0,0,0,0]),
           flag_entropy_instead_u::Int32=Int32(0),
           flag_doubleprecision::Int32=Int32(0),
           flag_ic_info::Int32=Int32(0),
           lpt_scalingfactor::Float32=Float32(0.),
           fill::Vector{Int32}=Int32.(zeros(12)))

          new(npart,
              massarr,
              time,
              z,
              flag_sfr,
              flag_feedback,
              nall,
              flag_cooling,
              num_files,
              boxsize,
              omega_0,
              omega_l,
              h0,
              flag_stellarage,
              flag_metals,
              npartTotalHighWord,
              flag_entropy_instead_u,
              flag_doubleprecision,
              flag_ic_info,
              lpt_scalingfactor,
              fill)
    end
end


"""
    head_to_obj(filename::String)

Returns the header of a snapshot as a `SnapshotHeader` object.
"""
function head_to_obj(filename)

    h = SnapshotHeader()

    f = open(filename)
    blocksize = read(f, Int32)

    if blocksize == 8
        swap = false
        snap_format = 2
    elseif blocksize == 256
        swap = false
        snap_format = 1
    else
        blocksize = bswap(blocksize)
        if blocksize == 8
            swap = true
            snap_format = 2
        elseif blocksize == 256
            swap = true
            snap_format = 1
        else
            error("incorrect file format encountered when reading header of $filename")
        end
    end

    #println("Reading snapshot format: ", snap_format)

    if snap_format == 2
        seek(f, 16)
        skip_line = read(f, Int32)
    end

    h.npart = read!(f, Array{Int32,1}(undef,6))
    h.massarr = read!(f, Array{Float64,1}(undef,6))
    h.time = read(f, Float64)
    h.z = read(f, Float64)
    h.flag_sfr = read(f, Int32)
    h.flag_feedback = read(f, Int32)
    h.nall = read!(f, Array{UInt32,1}(undef,6))
    h.flag_cooling = read(f, Int32)
    h.num_files = read(f, Int32)
    h.boxsize = read(f, Float64)
    h.omega_0 = read(f, Float64)
    h.omega_l = read(f, Float64)
    h.h0 = read(f, Float64)
    h.flag_stellarage = read(f, Int32)
    h.flag_metals = read(f, Int32)
    h.npartTotalHighWord = read!(f, Array{UInt32,1}(undef,6))
    h.flag_entropy_instead_u = read(f, Int32)
    h.flag_doubleprecision = read(f, Int32)
    h.flag_ic_info = read(f, Int32)
    h.lpt_scalingfactor = read(f, Float32)

    close(f)

    return h
end

"""
    head_to_dict(filename::String)

Returns the header of a snapshot as a dictionary.
"""
function head_to_dict(filename::String)

        header = Dict()

        f = open(filename)
        blocksize = read(f, Int32)

        if blocksize[1] == 8
            swap = 0
            snap_format = 2
        elseif blocksize[1] == 256
            swap = 0
            snap_format = 1
        else
            blocksize[1] = bswap(blocksize[1])
            if blocksize[1] == 8
                swap = 1
                snap_format = 2
            elseif blocksize[1] == 256
                swap = 1
                snap_format = 1
            else
                println("incorrect file format encountered when reading header of", filename)
            end
        end

        if snap_format == 2
            seek(f, 16)
            skip_line = read(f, Int32)
        end

        header["snap_format"] = snap_format
        header["PartTypes"] = ["PartType0", "PartType1", "PartType2",
                               "PartType3", "PartType4", "PartType5"]
        header["npart"] = read!(f, Array{Int32,1}(undef,6))
        header["massarr"] = read!(f, Array{Float64,1}(undef,6))
        header["time"] = read(f, Float64)
        header["redshift"] = read(f, Float64)
        header["flag_sfr"] = read(f, Int32)
        header["flag_feedback"] = read(f, Int32)
        header["nall"] = read!(f, Array{UInt32,1}(undef,6))
        header["flag_cooling"] = read(f, Int32)
        header["num_files"] = read(f, Int32)
        header["boxsize"] = read(f, Float64)
        header["omega_m"] = read(f, Float64)
        header["omega_l"] = read(f, Float64)
        header["hubble"] = read(f, Float64)
        header["flag_stellarage"] = read(f, Int32)
        header["flag_metals"] = read(f, Int32)
        header["npartTotalHighWord"] = read!(f, Array{UInt32,1}(undef,6))
        header["flag_entropy_instead_u"] = read(f, Int32)
        header["flag_doubleprecision"] = read(f, Int32)
        header["flag_ic_info"] = read(f, Int32)
        header["lpt_scalingfactor"] = read(f, Float32)

        close(f)

        return header

end


"""
    read_header(filename::String)

Reads the header of a snapshot file or file base (without .0, .1, etc.)
and returns a SnapshotHeader object.

See also: [`head_to_obj`](@ref)
"""
function read_header(filename::String)
    filename = select_file(filename, 0)
    return head_to_obj(filename)
end
