# GadgetIO.jl

This package provides IO functionality to work with the Smoothed-Particle Hydrodynamics code [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget/) by Volker Springel.

It is taylored for working with the development version of P-Gadget3, specifically OpenGadget3 developed by [Klaus Dolag](https://www.usm.uni-muenchen.de/~dolag/) and contributers. Development is focused on IO for Binary Format 2.
A lot of the routines are based on IDL scripts by Klaus Dolag (not public).

If you use `GadgetIO.jl` in your publications please cite [Böss & Valenzuela](https://zenodo.org/badge/latestdoi/270966661).

# Table of Contents

```@contents
Pages = [ "index.md",
          "install.md",
          "file_infos.md",
          "read_snapshots.md", 
          "read_subfind.md",
          "write_data.md",
          "api.md" 
        ]
Depth = 3
```
