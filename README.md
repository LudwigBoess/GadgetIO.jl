| **Documentation**                                                 | **Build Status**                                                                                | **License**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/GadgetIO.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/GadgetIO.jl/dev) | [![Build status](https://github.com/LudwigBoess/GadgetIO.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml/badge.svg)](https://github.com/LudwigBoess/GadgetIO.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml) [![codecov.io](https://codecov.io/gh/LudwigBoess/GadgetIO.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/GadgetIO.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) |


# GadgetIO.jl

This package is a subproject of [GadJet.jl](https://github.com/LudwigBoess/GadJet.jl) and provides some basic IO functionality to work with the Smoothed-Particle Hydrodynamics code [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget/) by Volker Springel.

It is taylored for working with the development version of P-Gadget3, specifically OpenGadget3 developed by Klaus Dolag. Development is focused on IO for Binary Format 2.

Please see the [Documentation](https://ludwigboess.github.io/GadgetIO.jl/dev/) for details.

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