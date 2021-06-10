Write Data
==========

```@meta
CurrentModule = GadgetIO
DocTestSetup = quote
    using GadgetIO
end
```

GadgetIO.jl can write snapshots that can be used as initial conditions.

Format 2
--------

The safest way to write snapshots is in Format 2.
Simply set up your header object and the arrays you want to write in the correct data format.
For the header this is the struct [`Header`](@ref) and for data its usually `Array{Float32,2}`.
You can then write an initial condition file by writing the header and the individual data blocks.

```julia
f = open(filename, "w")
write_header(f, header)
write_block( f,  pos, "POS")
write_block( f,  vel, "VEL")
write_block( f,  id,  "ID")
close(f)
```

Please note that you have to combine the arrays for individual particles in the correct order.

If you want to write the `INFO` block as well you can use

```julia
write_info_block(f, info)
```

where `info` is a `Vector` of [`InfoLine`](@ref).

Format 1
--------

Writing in format 1 works the same as above, but you need different function values.
Also you need to make sure the blocks are in the order gadget expects them to be!

```julia
f = open(filename, "w")
write_header(f, header, snap_format=1)
write_block( f,  pos,    snap_format=1)
write_block( f,  vel,    snap_format=1)
write_block( f,  id,     snap_format=1)
close(f)
```
