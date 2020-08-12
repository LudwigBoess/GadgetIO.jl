[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md)
[![Build Status](https://travis-ci.org/LudwigBoess/GadgetIO.jl.svg?branch=master)](https://travis-ci.org/LudwigBoess/GadgetIO.jl)
[![codecov.io](https://codecov.io/gh/LudwigBoess/GadgetIO.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/GadgetIO.jl?branch=master)


# GadgetIO.jl

This package is a subproject of [GadJet.jl](https://github.com/LudwigBoess/GadJet.jl) and provides some basic IO functionality to work with the SPH code "Gadget" by Volker Springel (doi:10.1111/j.1365-2966.2005.09655.x).

Documentation can be found [here](https://gadjetjl.readthedocs.io/en/latest/index.html).

Any help and contribution is greatly appreciated, as this is still a work in progress.

Quickstart
==========

Reading Data
------------

If you want to read a simulation snapshot into memory with GadJet.jl, it's as easy as this:

```julia
    data = read_snap(filename)
```

This will return a dictionary with the header information in `data["Header"]` and the blocks sorted by particle type.

As an example, this is how you would access the positions of the gas particles:

```julia
    data["Parttype0"]["POS"]
```

If you only want to read a specific block for a single particle type (e.g. positions of gas particles) you can use the function with a specified blockname and particle type like so:

```julia
    pos = read_snap(filename, "POS", 0)
```

This will return an array of the datatype of your simulation, usually Float32.