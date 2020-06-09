"""
    head_to_obj(filename::String)

Returns the header of a snapshot as a `Header` object.
"""
function head_to_obj(filename)

    h = Header()

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

Reads the header of a snapshot and returns a Header object.

See also: [`head_to_obj`](@ref)
"""
function read_header(filename::String)
    return head_to_obj(filename)
end
