"""
            Functions to write snapshots.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

function write_format2_block_header(f::IOStream, blockname::String)

    # fill blockname to 4 characters
    if length(blockname) < 4
        padding = 4 - length(blockname)
        for i=1:padding
            blockname = blockname * " "
        end
    end

    # transform the name into an array of charcters.
    # Since Julia uses 2 byte chars and C uses 1 byte I have to take
    # a different approach and use Int8s.
    write_name =  Int8.(collect(blockname))

    @info "Writing block: $blockname"

    # write blocksite (8 bytes)
    write(f, Int32(8))

    # write block name
    write(f, write_name)

    return f
end


"""
    write_block(f::IOStream, data,
                blockname::String="";
                snap_format::Integer=2)

Write `data` to a block to an opened file `f`.
"""
function write_block(f::IOStream, data,
                     blockname::String="";
                     snap_format::Integer=2)


    # turn array of arrays d into one large array
    #data = reduce(vcat, d)     # OLD! Keep for later, maybe

    # get total number of particles to write
    N = length(data[:,1])

    # write blocksize
    dtype = typeof(data[1,1])
    dims = length(data[1,:])
    blocksize = Int32(N * sizeof(dtype) * dims)

    if snap_format == 2

        # break if blockname not specified
        if blockname == ""
            error("Please specify blockname!")
        end

        # write header
        f = write_format2_block_header(f, blockname)

        # write end if name header, namely:
        # size to skip from directly after the block name and size of
        # the format 2 block header
        write(f, Int32(blocksize + Int32(8)))
        write(f, Int32(8))
    end

    # write blocksize
    write(f, blocksize)

    # write the block. Since Julia stores the arrays differently in memory
    # they have to be transposed before the can be written.
    write(f, copy(transpose(data)))

    @info "Writing block done."

    write(f, blocksize)
end


"""
    write_header(f::IOStream, h::Header; snap_format::Integer=2)

Writes the header block to an opened file `f`.
"""
function write_header(f::IOStream, h::Header; snap_format::Integer=2)

    if snap_format == 2
        f = write_format2_block_header(f, "HEAD")
        write(f, Int32(264))
        write(f, Int32(8))
    end

    # write blocksize
    write(f, Int32(256))

    # write header to file
    for fields in fieldnames(Header)

        write(f, getfield(h, fields))

    end

    # write blocksize
    write(f, Int32(256))
end


# Format 2 Blocknames
# "POS "
# "VEL "
# "ID  "
# "MASS"
# "U   "
# "RHO "
# "NE  "
# "NH  "
# "HII "
# "HeI "
# "HeII"
# "He3 "
# "H2I "
# "H2II"
# "HM  "
# "HD  "
# "DI  "
# "DII "
# "HeHp"
# "HSML"
# "VALP"
# "SFR "
# "AGE "
# "DETI"
# "HSMS"
# "ACRS"
# "Z   "
# "POT "
# "PHI "
# "PHID"
# "ACCE"
# "ENDT"
# "STRD"
# "STRO"
# "STRB"
# "SHCO"
# "TSTP"
# "BFLD"
# "BFSM"
# "DBDT"
# "VBLK"
# "VRMS"
# "VTAN"
# "VRAD"
# "TNGB"
# "DPP "
# "VDIC"
# "VROT"
# "VORT"
# "DIVB"
# "ACVC"
# "AMDC"
# "VTRB"
# "LTRB"
# "ADYN"
# "EDYN"
# "PHI "
# "XPHI"
# "GPHI"
# "ROTB"
# "SRTB"
# "EULA"
# "EULB"
# "COOR"
# "CONR"
# "DENN"
# "EGYP"
# "EGYC"
# "CRC0"
# "CRQ0"
# "CRP0"
# "CRE0"
# "CRn0"
# "CRco"
# "CRdi"
# "BHMA"
# "ACRB"
# "BHMD"
# "BHPC"
# "BHMB"
# "BHMI"
# "BHMR"
# "MACH"
# "SHSP"
# "SHRH"
# "SHPR"
# "SHVU"
# "SHCP"
# "SHNR"
# "DTEG"
# "PSCS"
# "PSDE"
# "PSEN"
# "PSXC"
# "DJMP"
# "EJMP"
# "CRDE"
# "TIPS"
# "DIPS"
# "CACO"
# "FLDE"
# "STDE"
# "SOMA"
# "PSDE"
# "ANRA"
# "LACA"
# "SHOR"
# "INDE"
# "aTEMP"  # remove a!
# "XNUC"
# "P   "
# "RADG"
# "RADA"
# "ET  "
# "PTSU"
# "DMNB"
# "NTSC"
# "SHSM"
# "SRHO"
# "SVEL"
# "DMHS"
# "DMDE"
# "DMVD"
# "VDMH"
# "VDMD"
# "Zs  "
# "CII "
# "CIa "
# "Cagb"
# "ZAge"
# "ZAlv"
# "iM  "
# "CLDX"
# "HOTT"
# "aTEMP"  # remove a!
# "TRCK"
# "ZsMT"
# "ZSMT"
# "DABN"
# "DABG"
# "DABS"
# "DCFB"
# "DCFA"
# "DCFD"
