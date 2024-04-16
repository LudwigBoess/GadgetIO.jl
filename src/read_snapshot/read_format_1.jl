"""
    read_block_format1(filebase, blockname::AbstractString; kwargs...)
    read_block_format1(filebase, blocknames::AbstractVector{<:AbstractString}; kwargs...)

Reads a block in a snapshot with given name. Defaults to reading gas particles. Block names are case sensitive.

# Keyword Arguments
- `parttype=0`: Which particle type to read (0-indexed)
    - `0`: Gas (default)
    - `1`: Dark Matter
    - `2`: Boundary
    - `3`: Bulge
    - `4`: Stars 
    - `5`: Black Holes
    - `-1`: All particle types
- `infolines::Union{Nothing,AbstractVector{InfoLine}}=nothing`: [`InfoLine`](@ref)s for the given file in the correct order. Needs to be used if the order of blocks is non-standard.
- `h::Union{Nothing,SnapshotHeader}=nothing`: `SnapshotHeader` for given file. Can be used to speed up IO.
- `is_inifile=false`: whether the file is an initial condition file or not
- `has_chem=false`: whether the simulation was run with chemistry (blocks `NE` and `NH`)
- `has_bh=false`: whether the simulation was run with black holes (blocks `BHMA` and `BHMD`)

# Default InfoLines
- `POS`: positions
- `VEL`: velocities
- `ID`: particle ID
- `MASS`: mass of particles
- `U`: internal energy of gas particles - if gas particles are present
- `RHO`: density - if gas particles are present
- `NE`: number density of free electrons - if `has_chem=true`
- `NH`: number density of neutral hydrogen - `has_chem=true`
- `HSML`: smoothing length of gas particles - if gas particles are present
- `SFR`: star formation rate - if `h.flag_sfr=1`
- `AGE`: Expansion factor or time at which star was born - if `h.flag_sfr=1`
- `BHMA`: true black hole mass (MASS contains the dynamical mass) - only if `has_bh=true`
- `BHMD`: black hole accretion rate - only if `has_bh=true`
- `DUMMY`: Unknown property for all particle types
"""
function read_block_format1(filebase, blockname::AbstractString; h::Union{Nothing,SnapshotHeader}=nothing, parttype=0, infolines::Union{Nothing,AbstractVector{InfoLine}}=nothing, is_inifile=false, has_chem=false, has_bh=false)
    if isnothing(h)
        h = read_header(filebase)
    end

    if isnothing(infolines)
        infolines = _get_default_infolines_format1(h, is_inifile, has_chem, has_bh)
    end

    infoind = findfirst(i -> i.block_name == blockname, infolines)
    @assert !isnothing(infoind) "Block $blockname is not known for format 1."

    infoline = infolines[infoind]

    n_to_read = _n_to_read_format1(h.nall, h.massarr, infoline, blockname, parttype)
    arr = allocate_data_array(infoline, n_to_read)

    next_ind = 1
    num_files = h.num_files

    for file in 0:num_files-1
        filename = select_file(filebase, file)

        htmp = file == 0 ? h : read_header(filename)
        n_to_read_tmp = _n_to_read_format1(htmp.npart, htmp.massarr, infoline, blockname, parttype)

        f = open(filename)

        # skip header
        seek(f, 256 + 2 * 4)

        for infoline in infolines
            if blockname == infoline.block_name
                # read array
                if infoline.n_dim > 1
                    _read_to_arr_format1!(f, @view(arr[1:infoline.n_dim, next_ind:next_ind+n_to_read_tmp-1]), infoline, parttype, htmp)
                else
                    _read_to_arr_format1!(f, @view(arr[next_ind:next_ind+n_to_read_tmp-1]), infoline, parttype, htmp)
                end

                break
            end

            _skip_block_format1(f, infoline, htmp)
        end

        close(f)

        next_ind += n_to_read_tmp
    end

    return arr
end

function read_block_format1(filebase, blocknames::AbstractVector{<:AbstractString}; h::Union{Nothing,SnapshotHeader}=nothing, parttype=0, infolines::Union{Nothing,AbstractVector{InfoLine}}=nothing, is_inifile=false, has_chem=false, has_bh=false)
    @assert 0 ≤ parttype ≤ 5 "A particle type has to be passed between 0 and 5"

    if isnothing(h)
        h = read_header(filebase)
    end

    if isnothing(infolines)
        infolines = _get_default_infolines_format1(h, is_inifile, has_chem, has_bh)
    end

    # prepare dictionary for particle storage
    d = Dict{String, VecOrMat{T} where T}()

    for blockname in blocknames
        d[blockname] = read_block_format1(filebase, blockname; h, parttype, infolines, is_inifile, has_chem, has_bh)
    end

    return d
end


function _n_to_read_format1(npart, massarr, infoline, blockname, parttype)
    is_present = infoline.is_present

    # all particles
    if parttype == -1
        if blockname == "MASS"
            return Int(sum(npart[massarr .== 0]))
        end

        return Int(sum(npart[is_present .== 1]))
    end

    # just a given particle type
    ind = parttype + 1

    if blockname == "MASS"
        @assert iszero(massarr[ind]) "The block MASS is not available for particle type $parttype. Please use `h.massarr[$(ind)]` instead (with `h = read_header(filename)`."
    end

    @assert is_present[ind] == 1 "The block $blockname is not available for particle type $parttype."

    return Int(npart[ind])
end


function _get_default_infolines_format1(h, is_inifile, has_chem, has_bh)
    @assert iszero(h.flag_doubleprecision) "GadgetIO is currently not able to read format 1 with double precision output. Please pass your own infolines as a keyword argument."

    d = default_info_lines
    infolines = get_info.((d,), ["POS", "VEL", "ID"])

    if any(iszero, h.massarr)
        push!(infolines, get_info(d, "MASS"))
    end

    if h.nall[1] > 0
        push!(infolines, get_info(d, "U"))

        if !is_inifile
            push!(infolines, get_info(d, "RHO"))

            if has_chem
                push!(infolines, get_info(d, "NE"))
                push!(infolines, get_info(d, "NH"))
            end

            push!(infolines, get_info(d, "HSML"))

            if h.flag_sfr == 1
                push!(infolines, get_info(d, "SFR"))
            end

            if h.flag_stellarage == 1
                push!(infolines, get_info(d, "AGE"))
            end

            if has_bh
                push!(infolines, get_info(d, "BHMA"))
                push!(infolines, get_info(d, "BHMD"))
            end
        end
    end

    return infolines
end

function _skip_block_format1(f, infoline, h)
    nbytes = read(f, UInt32)
    skip(f, nbytes)
    nbytes2 = read(f, UInt32)

    nbytessum = _n_to_read_format1(h.npart, h.massarr, infoline, infoline.block_name, -1) * infoline.n_dim * sizeof(infoline.data_type)

    @assert nbytes == nbytessum "It appears that the block $(infoline.block_name) is not what it was expected to be. It has $nbytes bytes instead of the expected $nbytessum."
    @assert nbytes == nbytes2 "Something is wrong with how the block $(infoline.block_name) is saved. The two byte sizes are $nbytes and $nbytes2 bytes."
end

function _read_to_arr_format1!(f, arr, infoline, parttype, h)
    # double check the correctness of this block
    nbytes = read(f, UInt32)
    nbytessum = _n_to_read_format1(h.npart, h.massarr, infoline, infoline.block_name, -1) * infoline.n_dim * sizeof(infoline.data_type)
    @assert nbytes == nbytessum "It appears that the block $(infoline.block_name) is not what it was expected to be. It has $nbytes bytes instead of the expected $nbytessum."

    # all particles
    if parttype == -1
        return read!(f, arr)
    end

    # just a given particle type
    ind = parttype + 1
    nskip = Int(sum(h.npart[findall(infoline.is_present[1:(ind-1)] .== 1)]))
    skip(f, nskip * infoline.n_dim * sizeof(infoline.data_type))

    return read!(f, arr)
end
