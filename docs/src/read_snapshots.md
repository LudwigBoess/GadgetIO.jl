# Read Snapshot Data


```@meta
CurrentModule = GadgetIO
DocTestSetup = quote
    using GadgetIO
end
```

---
**NOTE**

From v0.4 and up snapshots are read in proper column-major order, as it should be for Julia.
This means that position data for particle `i` neads be accessed as:
```
x = pos[1,i]
y = pos[2,i]
z = pos[3,i]
```
---


`GadgetIO.jl` is specialized to read `Gadget` snapshots of `Format 2`. The structure of a `Format 2` snapshot is as follows:

```
8              # size of the blockname block (Int32)
BLOCKNAME      # Blockname (4*Char)
8+SIZE_BLOCK   # number of bytes to skip if block should not be read
8              # end of blockname block

SIZE_BLOCK     # size of the current block in bytes
{...}          # content of the block ordered by particle type
SIZE_BLOCK     # end of the current block
```

which repeats for every block.

## Filename

Before we get started, a short bit for clarification: In what follows we need to distinguish between `filename` and `filebase`.
For small simulations the data is written to a single file. In that case you can simply supply an absolute, or relative path to the one snapshot file as `filename`.

For larger simulations the data may be distributed over multiple files, which again may be in individual snapshot directories. With the snapshots being distributed over multiple files you need to supply the base-name `filebase`. Assuming you want to read snapshot 140, which is in the snapshot directory 140 the `filebase` is

```julia
filebase = "path/to/your/snapshot/directories/snapdir_140/snap_140"
```

In the relevant functions `GadgetIO.jl` will then automatically loop through the sub-snapshots which end in ".0", ".1", ... , ".N".

## Reading a snapshot

There are multiple ways to read the data in the snapshot

## Full snapshot
If you want to read a simulation snapshot into memory with `GadgetIO.jl`, it's as easy as this:

```@docs
read_snap
```

This will return a dictionary with the header information in `data["Header"]` and the blocks sorted by particle type.

As an example, this is how you would access the positions of the gas particles:

```julia
    data["Parttype0"]["POS"]
```


## Specific blocks

If you only want to read a specific block for a single particle type, e.g. positions of gas particles, you can use the function `read_block` with a specified blockname and particle type.

```@docs
read_block
```

This will return an array of the datatype of your simulation, usually `Float32`.

If the snapshot has no info block this will fail unfortunately.

You can still read the specific block by supplying a hand-constructed [`InfoLine`](@ref) object:

```@docs
InfoLine
```

and passing that to the function [`read_block`](@ref):

```julia
pos = read_block(filename, "POS", info=pos_info, parttype=0)
```

where `pos_info` is an [`InfoLine`](@ref) object.


I will collect some example `InfoLine` objects in a later release to be able to read some common blocks even without an `INFO` block.

Since `v0.5` [`read_snap`](@ref) and [`read_block`](@ref) also work if you pass them a `file_base`.

Since `v0.7` [`read_block`](@ref) reads the full block if `parttype=-1` is set.

## Read Subvolumes

If you only want to read a subvolume of the whole simulation you can do this in two ways.
### Cubes

To get all particles within a cubic box of the simulation you can use the functions [`read_particles_in_box`](@ref) or [`read_particles_in_volume`](@ref).

[`read_particles_in_box`](@ref) takes a box defined by a lower-left corner and an upper-right corner and reads all requested blocks and particles in that volume.

```@docs
read_particles_in_box
```

You can define an array of blocks you want to read, these will be read into a dictionary.

[`read_particles_in_volume`](@ref) is a simple wrapper around [`read_particles_in_box`](@ref), where you can define a central position and a radius around it and it will construct the box containing that sphere for you and read all particles in it.

```@docs
read_particles_in_volume
```

### Arbitrary Geometries

For reading particles in more complex geometries you can use

```@docs
read_particles_in_geometry
```

You can use built-in geometries like [`GadgetCube`](@ref), [`GadgetSphere`](@ref) and [`GadgetCylinder`](@ref).

If you want to extend the functionality you can define your own geometry as

```julia
struct YourGeometry{T} <: AbstractGadgetGeometry
    prop1::Vector{T}
    prop2::Vector{T}
    prop3::T
    (...)
end
```

and define the functions [`get_geometry_box_corners`](@ref) and [`get_geometry_mask`](@ref).
[`get_geometry_box_corners`](@ref) has to return a `Tuple` of two vectors which define the lower left and upper right corner of a box that contains the `geometry`. [`get_geometry_mask`](@ref) has to return an array of indices for which `pos` is contained in the `geometry`.

In all functions `parttype` defines the particle type to be read, as in the previous read functions and `verbose` gives console output.
There are also multiple dispatch versions of all functions available that only take a single `block` as input and return an array with the values instead of a dictionary.



## Peano-Hilbert key based reading

For large simulations Gadget distributes snapshots over multiple files. These files contain particles associated with specific Peano-Hilbert keys.

If you call [`read_particles_in_box`](@ref) or [`read_particles_in_volume`](@ref) with the keyword argument `use_keys=true` (which is the default case) it constructs the peano hilbert keys, selects the relevant files and reads the particles from these files into a dictionary. This is considerably faster than the brute-force attempt.


## Brute-Force Reading
If you call [`read_particles_in_box`](@ref) or [`read_particles_in_volume`](@ref) with the keyword argument `use_keys=false` it reads all particles over all distributed files which are contained in the requested subvolume.
This takes quite a lot longer than the key based reading, but sometimes it's the only option.
To speed this up you can apply the filtering only once and store the [Read positions](@ref).



## Custom Filtering

If you want to read multiple blocks in a simulation whose snapshots have been distributed over a number of sub-snapshots you can use [`read_blocks_filtered`](@ref).

```@docs
read_blocks_filtered
```
This will read the specified `blocks` for all particles that pass the `filter_function`. This can be useful if you don't know where the region you are interested in is located and don't have enough memory to read in all particles.
Alternatively you can also provide a `Dict` with `read_positions` as described in [Read positions](@ref).

### Filter functions

The `filter_function` can be any function that takes a `String` input and returns an `Array` of `Integer`s, or `CartesianCoordinates`.
For example, if you want to filter all particles with a Mach number larger than 1:

```julia
function mach_gt_1(snap_file)
    mach = read_snap(snap_file, "MACH", 0)
    sel  = findall( mach .> 1.0 )
    return sel
end
```

Or if you want to trick the function into reading all particles after all:

```julia
function pass_all(snap_file)
    h = read_header(snap_file)
    return collect(1:h.npart[1])
end
```

As an example, to read positions, velocity and ID of all shocked particles from distributed snapshots use

```julia
blocks = ["POS", "VEL", "ID"]
data = read_blocks_filtered(snap_base, blocks, filter_function=mach_gt_1, parttype=0)
```

Just as a reminder from above you can read single blocks into an array by using [`read_snap`](@ref) and [`read_block`](@ref):

```julia
pos = read_block(snap_base, "POS", parttype=0)
```

This requiers `GadgetIO.jl` v0.5 though.

### Read positions

To avoid having to filter all files each time you want to read a snapshot you can also split the steps.
You can first filter the particles to find the positions of the particles within the data blocks with

```@docs
find_read_positions
```
and then save the result as a binary file with

```@docs
save_read_positions
```

To re-use the `read_positions` you can load them from file using [`load_read_positions`](@ref)

```julia
read_positions = load_read_positions(save_file)
```

This can then be used to read any number of blocks with

```julia
blocks = ["POS", "VEL", "ID"]
data = read_blocks_filtered(snap_base, blocks, read_positions=read_positions, parttype=0)
```

## Reading particles by referenced ID

If you want to select specific particles to read from an array of `IDs` you can do this with

```@docs
read_particles_by_id
```

`snap_base` defines the target snapshot, or the snapshot basename, `selected_ids` contains the list of IDs of the particles you want to read and `blocks` containes the blocknames of the blocks you want to read.
If the simulation is too large to read the whole snapshot into memory you can give values for `pos0` and `r0` to read only a specific region with [`read_particles_in_volume`](@ref). See [Read Subvolumes](@ref) for details on this.

This will return a dictionary with all requested blocks.



## Example

If you want to, e.g. read positions, velocities, masses, density and hsml for all gas particles within the virial radius of the most massive halo of a simulation you can do this as follows.

Assuming `pos_halo` is the position of the center of mass of the halo and `r_vir` is its virial radius you read the data with

```julia
blocks = ["POS", "VEL", "MASS", "RHO", "HSML"]

data   = read_particles_in_volume(filename, blocks, pos_halo, r_vir,
                                  parttype=0,
                                  verbose=true)
```

This will return a dictionary with the blocks as keys and containing the arrays for the particles.

```julia
data["POS"]  # array of positions
data["RHO"]  # array of densities
(...)
```