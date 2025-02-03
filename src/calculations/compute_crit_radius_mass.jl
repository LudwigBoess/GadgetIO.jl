

function get_radius(pos::Matrix{<:Real}, pos0::Vector{<:Real})

    # allocate memory for the radius
    r = Vector{Float64}(undef, size(pos, 2))

    # loop over all particles
    for i in 1:size(pos, 2)
        r0 = 0.0
        for dim = 1:3
            r0 += (pos[dim, i] - pos0[dim])^2
        end
        r[i] = sqrt(r0)
    end

    return r
end

function compute_Rcrit_Mcrit(halo::HaloID, snap_base::String, sub_base::String, rho_crit::Real;
                            radius::Union{Real,Nothing}=nothing,
                            rad_scale::Real=10.0, halo_type::Integer=1,
                            parttype::Integer=0, verbose::Bool=true,
                            use_keys::Bool=true)
    
    # select subfind file to read
    sub_file = select_file(sub_base, halo.file)
    h = read_header(sub_file)

    if verbose
        @info "Reading IDs in halo..."
        flush(stdout)
        flush(stderr)
        t1 = Dates.now()
    end

    # to compute Rcrit and Mcrit we need positions and masses
    blocks = ["POS", "MASS"]

    # read all IDs of the particles contained in a halo
    data = read_particles_in_halo(snap_base, blocks, sub_base, halo; 
                                  radius, rad_scale, halo_type, parttype, 
                                  verbose, use_keys)

    if verbose
        t2 = Dates.now()
        @info "Data read. Took: $(t2 - t1)"
        println()
        @info "Done!"

        @info "Computing radial profile..."
        flush(stdout)
        flush(stderr)
        t1 = Dates.now()
    end

    # position of halo
    pos_block = get_pos_block_name(halo_type)
    pos0      = read_subfind(sub_file, pos_block)[:,halo.id]

    r = get_radius(data["POS"], pos0)
    sorted = sortperm(r)

    # store sorted arrays for radius and mass
    r = r[sorted]
    M = data["MASS"][sorted]

    # compute the radial Mass profile
    Msum = Vector{Float64}(undef, length(data["MASS"]))
    Msum[1] = M[1] 
    for i in 2:length(r)
        Msum[i] = Msum[i-1] + M[i]
    end

    if verbose
        t2 = Dates.now()
        @info "Data read. Took: $(t2 - t1)"
        println()
        @info "Done!"

        @info "Computing critical point..."
        flush(stdout)
        flush(stderr)
        t1 = Dates.now()
    end

    # compute the critical radius
    crit_point = findfirst(Msum ./ (4pi/3 .* r.^3) .< rho_crit)
    Rcrit = r[crit_point]
    Mcrit = Msum[crit_point]

    if verbose
        t2 = Dates.now()
        @info "Done. Took: $(t2 - t1)"
        println()
        @info "Done!"
        flush(stdout)
        flush(stderr)
    end

    return Rcrit, Mcrit

end
