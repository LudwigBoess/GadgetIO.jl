"""
    get_total_particles(h::AbstractGadgetHeader, parttype::Integer)

Calculates to total number of particles present in the simulation. Accounts for integer overflow.
"""
function get_total_particles(npart::Integer, high_bit::Integer)
    high_bit * 2^32 + npart
end


"""
    get_total_particles(h::AbstractGadgetHeader, parttype::Integer)

Calculates to total number of particles present in the simulation. Accounts for integer overflow.
"""
function get_total_particles(h::SnapshotHeader, parttype::Integer)
    get_total_particles(h.nall[parttype+1], h.npartTotalHighWord[parttype+1])
end