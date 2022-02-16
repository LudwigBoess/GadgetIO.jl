# Read Subfind Data

`Gadget` contains an on-the-fly halo-finder as described in [Springel et al (2001)](https://ui.adsabs.harvard.edu/link_gateway/2001MNRAS.328..726S/doi:10.1046/j.1365-8711.2001.04912.x) or [Dolag et al (2009)](https://ui.adsabs.harvard.edu/link_gateway/2009MNRAS.399..497D/doi:10.1111/j.1365-2966.2009.15034.x).
This sections provides an overview of the functions you can use to work with this output.
Please note that you need to compile `Gadget` with `WRITE_SUB_IN_SNAP_FORMAT` to use this functionality.

## Reading the header

As in the normal snapshot the subfind output also contains a `HEAD` block with useful information.
You can read the header of the subfind output into a [`SubfindHeader`](@ref) object

```julia
struct SubfindHeader
    nhalos::Int32                       # number of halos in the output file
    nsubhalos::Int32                    # number of subhalos in the output file
    nfof::Int32                         # number of particles in the FoF
    ngroups::Int32                      # number of large groups in the output file
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    tothalos::UInt32                    # total number of halos over all output files
    totsubhalos::UInt32                 # total number of subhalos over all output files
    totfof::UInt32                      # total number of particles in the FoF
    totgroups::UInt32                   # total number of large groups over all output files
    num_colors::Int32                   # number of colors
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32
end
```

using

```julia
h = read_subfind_header(filename::String)
```

## Reading the subfind files


For convenience you can use a helper function provided by `GadgetIO` to read the block of the subfind output. Since each of the blocks is only relevant for either halos, subhalos, Fof or large groups you don't need to define a particly type, aka halo type in this case.

So in order to read the virial radius of the halos in a file you can simply use

```julia
R_vir = read_subfind(filename, "RVIR")
```

## Filtered read-in

If you want to read specific halos from the subfind output you can use the function

```julia
filter_subfind(filebase::String, blockname::String, filter_function::Function, nfiles::Integer=1)
```

This will return an array of [HaloID](@ref)s 

```julia
struct HaloID
    file::Int64
    id::Int64
end
```

You can use this to read specific files with [`read_subfind`](@ref) and select the array entry via the `id` field.

The `filter_function` argument takes any function that takes one input argument and returns `true` if the requirement is fulfilled, or `false` if not.

So to find e.g. all halos with a virial mass larger than ``10^{15} M_\odot`` you can use

```julia
find_mass_gt_1e15(M) = ( (M > 1.e15) ? true : false )

filtered_subfind = filter_subfind(filebase, "MVIR", find_mass_gt_1e15)
```

This will search all subfind files in the directory and check if they fulfill the filter. Since halos in subfind files are sorted by mass you can also supply a number of files to search with the argument `nfiles`. That way, if you are looking for a very massive halo you can constrian the reading and filtering to only the first `N` files to save time.

If you want more complex filtering you can also just pass a `filter_function` that takes a filename as an input and returns an `Array` of `Integer`s, or `CartesianCoordinates`, like was the case in [Filter functions](@ref).

```julia
filtered_subfind = filter_subfind(filebase, filter_function)
```

## Reading halo property by HaloID

You can read any property of the halo that passed the `filter_function` (see [Filtered read-in](@ref)) by using [`read_halo_prop`](@ref). So if you want to read e.g. the virial radius for the first halo that passed your `filter_function`

```julia
mvir = read_halo_prop(filebase, filtered_subfind[1], "MVIR")
```

HaloIDs can also be obtained in a vector of all halos by setting the keyword parameter `return_haloid` to `true`:

```julia
mvir, haloids = read_subfind(filename, "MVIR"; return_haloid=true)
```

## Reading halo properties by global halo index

If you have a global index of a halo from subfind (0-indexed, increasing over all subfiles of the subfind outputs), you can also read a halo's properties from subfind. To read a halo property from such a global halo index and convert it to a `HaloID` you can use [`read_halo_prop_and_id`](@ref). To only obtain the property, [`read_halo_prop`](@ref) can be used.
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


## Reading particles in a halo

If you want to read all particles associated with a FoF halo you can do this with the function [`read_particles_in_halo`](@ref)

```julia
read_particles_in_halo( snap_base::String,   blocks::Array{String},
                        sub_base::String,    halo::HaloID; 
                        rad_scale::Real=1.0, halo_type::Integer=1,
                        parttype::Integer=0, verbose::Bool=true)
```

This reads all blocks defined in `blocks` for the `halo` into a dictionary. `snap_base` and `sub_base` should point to the snap and subfind filebase as in other functions, or the files if you only have one file. The `HaloID` point to the selected halo. `halo_type` should be set to `1` for halos and `2` for subhalos. 
If you didn't get all the particles you were looking for it might be that the search radius for the read-in was too small. `rad_scale` defines the multiplication factor for the search radius. For halos the default search radius is ``r_{200}`` and for subhalos it's the half-mass radius.

