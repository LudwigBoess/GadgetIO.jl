# Finding the Main Halo

To work with a zoom simulation you usually need the centre of its main halo, for example to read or map the particles around it.
If the simulation ran with `Subfind` you can use [`find_most_massive_halo`](@ref) on its output, see [Read Subfind Data](@ref).
Without halo finder output you can find the main halo directly in the snapshot with

```@docs
find_main_halo
```

## Why not the potential minimum

The particle with the lowest potential in a zoom simulation is not necessarily in the main halo.
Gas particles that left the high-resolution region, and high-resolution particles next to massive low-resolution particles, can sit in a deeper potential than the centre of the main halo.
`find_main_halo` therefore first identifies the main halo from the high-resolution dark matter alone and only then looks for the potential minimum among the particles of that halo.

## How it works

1. A friends-of-friends (FoF) algorithm groups the high-resolution dark matter particles (particle type 1).
   By default the linking length is `b = 0.2` times the mean particle separation at the mean dark matter density, as in the FoF of OpenGadget3.
   The mean matter density is the total mass of all particles in the box divided by the box volume, so it does not depend on the unit system.
   Dark matter takes the fraction ``(\Omega_0 - \Omega_b) / \Omega_0`` of it.
   Since ``\Omega_b`` is not stored in the snapshot header, pass it with `omega_b` if it differs from the default of `0.04`.
2. The main halo is the group with the largest dark matter mass.
3. The centre of the main halo is the position of its particle with the lowest potential (`center=:potential`), which needs the `POT` block in the snapshot.
   For snapshots without potential use `center=:shrinking_sphere`, the shrinking-sphere centre of [Power et al. (2003)](https://ui.adsabs.harvard.edu/abs/2003MNRAS.338...14P) of the particles in the halo.

Before linking, the positions are unwrapped around one of the high-resolution particles, so a main halo that lies across the border of the periodic box stays in one piece.
The returned centre is then the periodic image next to the high-resolution region.

```julia
snap_base = "path/to/your/snapshot/directories/snapdir_140/snap_140"

# main halo with the potential minimum as centre
center = find_main_halo(snap_base)

# read the gas in a cube of 4 Mpc/h around it
data = read_particles_in_geometry(snap_base, ["POS", "MASS"], GadgetCube(center, 4000.0), parttype=0)
```

## Comparison to the on-the-fly halo finder

With the same linking length the FoF groups contain the same dark matter particles as the FoF groups of the halo finder in `Gadget`, which in addition attaches particles of the other types to them.
For the test snapshot of this package all 17 groups in the `Subfind` output are reproduced particle by particle.
Older simulations may use a different linking length than the current OpenGadget3 definition, the test snapshot for example 42 instead of 49.6 in code units.
In that case pass the linking length of the simulation with `linking_length` to reproduce its groups.

The main halo of `find_main_halo` is the FoF group with the largest dark matter mass, while [`find_most_massive_halo`](@ref) selects the halo with the largest spherical overdensity mass in the `Subfind` output.
At high redshift, when the progenitors of the main halo are still separate, the two can select different halos.

## Performance

The FoF sorts the particles into cells that are smaller than the linking length, so all particles in a cell are linked, and only neighbouring cells have to be searched for links.
Its run time grows about linearly with the number of particles: the 1.6×10⁷ high-resolution dark matter particles of a zoom simulation take about 7 s on a single core.
The cells are indexed by a dense two-dimensional array, which is small for the compact high-resolution region of a zoom simulation, but not for the full box of a cosmological simulation.

You can also use the FoF and the shrinking-sphere centre on their own:

```@docs
GadgetIO.fof_labels
GadgetIO.shrinking_sphere_center
GadgetIO.fof_linking_length
```
