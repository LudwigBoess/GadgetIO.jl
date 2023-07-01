| **Documentation**                                                 | **Build Status**                                                                                | **License**                                                                                | **Citation**
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/GadgetIO.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/GadgetIO.jl/dev) | [![Build status](https://github.com/LudwigBoess/GadgetIO.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml/badge.svg)](https://github.com/LudwigBoess/GadgetIO.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml) [![codecov.io](https://codecov.io/gh/LudwigBoess/GadgetIO.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/GadgetIO.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) | [![DOI](https://zenodo.org/badge/270966661.svg)](https://zenodo.org/badge/latestdoi/270966661) |


# GadgetIO.jl

This package provides some basic IO functionality to work with the Smoothed-Particle Hydrodynamics code [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget/) by Volker Springel.

You can use it to read/write snapshot and subfind output and we also provide some functionality to read the log files.

It is taylored for working with the development version of P-Gadget3, specifically OpenGadget3 developed by [Klaus Dolag](https://www.usm.uni-muenchen.de/~dolag/) and contributers. Development is focused on IO for Binary Format 2.
A lot of the routines are based on IDL scripts by Klaus Dolag (not public).

If you use `GadgetIO.jl` in your publications please cite [Böss & Valenzuela](https://zenodo.org/badge/latestdoi/270966661).

Please see the [Documentation](https://ludwigboess.github.io/GadgetIO.jl/dev/) for details.

Quickstart
==========

Reading Snapshot Data
------------

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

Gadget defines 6 particle types:
- `0`: Gas particles
- `1`: DM particles
- `2`: Disk/boundary particles
- `3`: Bulge/boundary particles
- `4`: Star particles
- `5`: Black Hole particles

So if you want to read for example the poistions of gas particles you can do this by using:

julia```
filename = "path/to/your/snapshot"

gas_pos = read_block(filename, "POS", parttype=0)
```

Similar for DM particles:
julia```
filename = "path/to/your/snapshot"

dm_pos = read_block(filename, "POS", parttype=1)
```

We recommend running Gadget with the compiler flag `WRITE_INFO_BLOCK` to make IO easier, however you can also read most blocks out of the box due to fall-back `InfoLine`s.
If you work on a development version of `Gadget` that contains output blocks that don't have fall-back `InfoLine`s you can also supply your own to the `read_block` function.

You can also read subvolumes of simulations, either along the Peano-Hilbert curve or in a brute force way.
Or you can read particles by custom filtering and store the memory-mapping to disk.
See the [Documentation](https://ludwigboess.github.io/GadgetIO.jl/dev/) for details.