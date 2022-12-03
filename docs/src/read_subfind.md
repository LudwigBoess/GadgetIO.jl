# Read Subfind Data

`Gadget` contains an on-the-fly halo-finder as described in [Springel et al (2001)](https://ui.adsabs.harvard.edu/link_gateway/2001MNRAS.328..726S/doi:10.1046/j.1365-8711.2001.04912.x) or [Dolag et al (2009)](https://ui.adsabs.harvard.edu/link_gateway/2009MNRAS.399..497D/doi:10.1111/j.1365-2966.2009.15034.x).
This sections provides an overview of the functions you can use to work with this output.
Please note that you need to compile `Gadget` with `WRITE_SUB_IN_SNAP_FORMAT` to use this functionality.

Please note that since subfind files are equivalent to snapshot files you can use the same functions as in if you [Read Snapshot Data](@ref).

## Reading the header

As in the normal snapshot the subfind output also contains a `HEAD` block with useful information.
You can read the header of the subfind output into a [`SubfindHeader`](@ref) object

```@docs
SubfindHeader
```

by using

```@docs
read_subfind_header
```

## Reading the subfind files

For convenience you can use a helper function provided by `GadgetIO.jl` to read the block of the subfind output. Since each of the blocks is only relevant for either halos, subhalos, FoF or large groups you don't need to define a particly type, aka halo type in this case.

```@docs
read_subfind
```

## Filtered read-in

If you want to read specific halos from the subfind output you can use the function

```@docs
filter_subfind
```

This will return an array of [HaloID](@ref)s 

```@docs
HaloID
```

This can be used either with [Reading halo property by HaloID](@ref) 

## Saving/Loading HaloIDs

If you do complex filtering and want to save the result you can use

```@docs
save_halo_ids
```

and load them the next time with 

```@docs
load_halo_ids
```


If you want to use [`HaloID`](@ref)s with [`read_block_filtered`](@ref) you can convert convert the [`HaloID`](@ref)s to `read_positions` by using

```@docs
halo_ids_to_read_positions
```


## Reading halo property by HaloID

You can read any property of the halo that passed the `filter_function` (see [Filtered read-in](@ref)) by using 

```@docs
read_halo_prop(::Any, ::AbstractString, ::HaloID)
```

So if you want to read e.g. the virial radius for the first halo that passed your `filter_function`

```julia
mvir = read_halo_prop(filebase, filtered_subfind[1], "MVIR")
```

HaloIDs can also be obtained in a vector of all halos by setting the keyword parameter `return_haloid` to `true`:

```julia
mvir, haloids = read_subfind(filename, "MVIR"; return_haloid=true)
```

## Reading halo properties by global halo index

If you have a global index of a halo from subfind (0-indexed, increasing over all subfiles of the subfind outputs), you can also read a halo's properties from subfind. To read a halo property from such a global halo index and convert it to a `HaloID` you can use 

```@docs
read_halo_prop_and_id
```

To only obtain the property use

```@docs
read_halo_prop(::Any, ::AbstractString, ::Integer)
```


So for example if you have the global halo id `i_global` and want to read the corresponding virial mass you can use one of the two following lines (note that it is faster to read properties via the HaloID since only a single file has to be read for that):

```julia
mvir, halo_id = read_halo_prop_and_id(filebase, i_global, "MVIR")
mvir = read_halo_prop(filebase, i_global, "MVIR")
```

To read the properties of multiple halos for which the halo indices are available, use one of the following two lines (to read the virial masses of the first four halos in subfind):

```julia
mvir = read_subfind(filebase, "MVIR", [0, 1, 2, 3])
mvir, haloids = read_subfind(filebase, "MVIR", [0, 1, 2, 3]; return_haloid=true)
```

The results are returned in the order of the given indices.

## Converting global halo indices to HaloIDs

If you want to convert global halo indices to [`HaloID`](@ref)s use

```@docs
global_idxs_to_halo_id
```

You can then save the `HaloID`s with [Saving/Loading HaloIDs](@ref) or convert them to `read_positions` with [`halo_ids_to_read_positions`](@ref).


## Reading particles in a halo

If you want to read all particles associated with a FoF halo you can do this with the function 

```@docs
read_particles_in_halo
```

This reads all blocks defined in `blocks` for the `halo` into a dictionary. `snap_base` and `sub_base` should point to the snap and subfind filebase as in other functions, or the files if you only have one file. The `HaloID` point to the selected halo. `halo_type` should be set to `1` for halos and `2` for subhalos. 
If you didn't get all the particles you were looking for it might be that the search radius for the read-in was too small. `rad_scale` defines the multiplication factor for the search radius. For halos the default search radius is ``r_{200}`` and for subhalos it's the half-mass radius.

