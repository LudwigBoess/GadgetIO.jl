# File information

```@meta
CurrentModule = GadgetIO
DocTestSetup = quote
    using GadgetIO
end
```

Since `Gadget` outputs in binary format it can be quite tedious to see what is actually contained in the file. For this `GadgetIO` provides a number of helper functions.

## Reading the header

The header block (`HEAD`) contains a number of informations about the file, for example how many particles of which kind are contained in the file, the current time or redshift, the cosmological parameters and so on.

The most convenient way to deal with this data is to read the header block into a [`SnapshotHeader`](@ref) `struct` by using

```@docs
read_header
```

the struct contains the following fields:

```@docs
SnapshotHeader
```

[`read_header`](@ref) is a wrapper around:

```@docs
head_to_struct
```

If you want to read the header information into a dictionary you can use:

```@docs
head_to_dict
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
get_block_positions
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

```@docs
InfoLine
```

You can use this as well to speed up read-in by passing the relevant `Array` entry to the keyword argument `info` in [`read_block`](@ref):

```julia
gas_pos = read_block(filename, "POS", parttype=0, info=info[1])
```

This can of course be combined with reading the block positions for additional speedup.

If there is no `INFO` block in the simulation you need to construct the `InfoLine` struct yourself and pass it to [`read_block`](@ref) like in the previous example.

