function assign_mass_from_header!(block, filename, h, parttype)
    if parttype != -1
        block .= h.massarr[parttype+1]
    else
        # assign masses for all particles
        n_read = 0
        # assign mass from header one by one
        for ptype = 0:5
            n_to_read = h.npart[ptype+1]
            if !iszero(n_to_read)
                if !iszero(h.massarr[ptype+1])
                    block[n_read+1:n_read+n_to_read] .= h.massarr[ptype+1]
                else
                    block[n_read+1:n_read+n_to_read] = read_block(filename, "MASS", parttype=ptype)
                end 
                n_read += n_to_read
            end
        end
    end

    return block
end