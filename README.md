# Quaycle-Elastic.jl

- This is a trimmed version of [Quaycle.jl](https://github.com/shipengcheng1230/Quaycle.jl) with only elastic part. Full building, testing and docs building are not included here.

- Documents & Examples: [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://shipengcheng1230.github.io/Quaycle.jl/dev)

- Status: ![Not Maintained](https://img.shields.io/badge/Repo%20Status-Achieved-yellow)

- Building: [![Build Status](https://travis-ci.com/shipengcheng1230/Quaycle-Elastic.jl.svg?branch=master)](https://travis-ci.com/shipengcheng1230/Quaycle-Elastic.jl)

- Citation: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3836243.svg)](https://doi.org/10.5281/zenodo.3836243)



# JGR2020

This repo constains the sample scripts for reproducing simulations and the outputs of catalogues.

- scripts:
    - `s01-domain.jl`: create a fault domain
    - `s02-greensfunc.jl`: compute dislocation-stress Green's function
    - `s03-parameters.jl`: set fault parameters
    - `s04-solve.jl`: solve models
    - `scanfunc.jl`: functions for computing event catalogues from raw output data
    - `solve.job`: an example sbatch script for solving each model using one node (multi-threading within node, array job)


- outputs (for catalogue only):

    All model outputs are in HDF5. `ic{#}/` denotes different initial conditions. `extend/` contains extended simulations (up to 1200 or 1800 years, end-year is appened in those output names).

    Catalog output names are like `resf-{parameter group}{value}.h5`. LF denote left half fault, RF right half fault. Fields in the catalogue outputs are:

    - `t`: time step in second
    - `ti`: time step in second (downsampled if `_stride != 1`)
    - `maxva`: max velocity in m/s (LF)
    - `maxvb`: max velocity in m/s (RF)
    - `maxvf`: max velocity in m/s (whole fault)
    - `mwa`: moment magnitude (LF)
    - `mwb`: moment magnitude (RF)
    - `ixba`: event start time (LF), correspond to `ti`
    - `ixbb`: event start time (RF), correspond to `ti`
    - `ixea`: event end index (LF), correspond to `ti`
    - `ixeb`: event end index (RF), correspond to `ti`
    - `riax`: event initial rutpure along strike position (m) (LF)
    - `ribx`: event initial rutpure along strike position (m) (RF)
    - `riaz`: event initial rutpure down dip position (m) (LF)
    - `ribz`: event initial rutpure down dip strike position (m) (RF)


    Slip ratio output names are like `srf-{parameter group}{value}.h5`. Fields in the catalogue outputs are:

    - `sr`: seismic ratio
    - `ar`: afterslip ratio


    `{output type}-s-{parameter group}{value}.h5` denote single VW model.

    ODEs solution output names are like `otf-{parameter group}{value}.h5` (not included in this repo).
    They contains three fields:

    - `t`: time step (s)
    - `v`: velocity (m/s)
    - `θ`: state variable (s)
    - `δ`: displacement (m)