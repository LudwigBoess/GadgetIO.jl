# File information

Since `Gadget` outputs in binary format it can be quite tedious to see what is actually contained in the file. For this `GadgetIO` provides a number of helper functions.

## Reading the header

The header block (`HEAD`) contains a number of informations about the file, for example how many particles of which kind are contained in the file, the current time or redshift, the cosmological parameters and so on.

The most convenient way to deal with this data is to read the header block into a [`SnapshotHeader`](@ref) `struct` by using

```@docs
read_header
```

the struct contains the following fields:

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
| `npartTotalHighWord::Vector{UInt32}` | If Npart > 1584^3 (>2^32) this contains a high bit: ntotal = npartTotalHighWord*2^32 + nall  |
| `flag_entropy_instead_u::Int32`      | 1 if snapshot U field contains entropy instead of internal energy, else 0              |
| `flag_doubleprecision::Int32`        | 1 if snapshot is in double precision, else 0                                           |
| `flag_ic_info::Int32`                | 1 if initial snapshot file contains an info block, else 0                              |
| `lpt_scalingfactor::Float32`         | factor to use second order ic generation                                               |
| `fill::Vector{Int32}`                | the HEAD block needs to be filled with zeros to have a size of 256 bytes               |


[`read_header`](@ref) is a wrapper around:

```@docs
GadgetIO.head_to_struct
```

If you want to read the header information into a dictionary you can use:

```@docs
GadgetIO.head_to_dict
```

## Getting the block names


```@docs
print_blocks
```


## Checking if a block is present


```@docs
block_present
```

## Getting the block positions within a file

If you want to read in multiple blocks into individual variables you can speed the process up significantly by first getting the starting positions of the blocks with

```@docs
GadgetIO.get_block_positions
```

This returns a dictionary with the block names as keys and the starting position in bytes as an `Integer`:

```julia
julia> block_pos["POS"]
1234567
```

and can be used with [`read_block`](@ref) by passing it to the keyword argument `block_position`:

```julia
gas_pos = read_block(filename, "POS", parttype=0, block_position=block_pos["POS"])
gas_vel = read_block(filename, "VEL", parttype=0, block_position=block_pos["VEL"])
gas_rho = read_block(filename, "RHO", parttype=0, block_position=block_pos["RHO"])
```

Please note that this only makes sense if you plan to read in multiple blocks. If the keyword argument is left out [`read_block`](@ref) will look for the block positions by itself.

## Reading the INFO block

If you compiled `Gadget` with `WRITE_INFO_BLOCK` the snapshot contains a block `INFO` that holds information on the name, datatype, dimensionality and presence per particle for each block. This simplifies read-in greatly and is always recommended!
You can obtain this block by using:

```@docs
read_info
```

This returns an `Array{InfoLine,1}`, where [`InfoLine`](@ref) is:

```
struct InfoLine([  block_name="", data_type=Float32, n_dim=Int32(0),
                is_present=zeros(Int32, 6) ])
```

Contains the data of a single entry in the `INFO` block of a Gadget snapshot.
The data can be accessed via the fields:

| Name                                 | Meaning                                                                                |
| :----------------------------------  | :------------------------------------------------------------------------------------- |
| `block_name::String`                 | name of the data block, e.g. "POS"                                                     |
| `data_type::DataType`                | datatype of the block, e.g. Float32 for single precision, Float64 for double           |
| `n_dim::Int32`                       | number of dimensions of the block, usually 1 or 3                                      |
| `is_present::Vector{Int32}`          | array of flags for which particle type this block is present,                          |
|                                      |  e.g. gas only:  [ 1, 0, 0, 0, 0, 0 ],                                                 |
|                                      |  or gas + BHs: [ 1, 0, 0, 0, 0, 1 ]                                                    |


You can use this as well to speed up read-in by passing the relevant `Array` entry to the keyword argument `info` in [`read_block`](@ref):

```julia
gas_pos = read_block(filename, "POS", parttype=0, info=info[1])
```

This can of course be combined with reading the block positions for additional speedup.

If there is no `INFO` block in the simulation you need to construct the `InfoLine` struct yourself and pass it to [`read_block`](@ref) like in the previous example.

