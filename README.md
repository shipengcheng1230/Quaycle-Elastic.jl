# Quaycle-Elastic.jl

- This is a trimmed version of [Quaycle.jl](https://github.com/shipengcheng1230/Quaycle.jl) with only elastic part. Full building, testing and docs building are not included here.

- Documents & Examples: [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://shipengcheng1230.github.io/Quaycle.jl/dev)

- Status: ![Not Maintained](https://img.shields.io/badge/Repo%20Status-Achieved-yellow)

- Building: [![Build Status](https://travis-ci.com/shipengcheng1230/Quaycle-ElasticPart.jl.svg?branch=master)](https://travis-ci.com/shipengcheng1230/Quaycle-ElasticPart.jl)


# JGR2020

This repo constains the sample scripts for reproducing simulations and the outputs of catalogues.

- scripts:
    - `s01-domain.jl`: create a fault domain
    - `s02-greensfunc.jl`: compute dislocation-stress Green's function
    - `s03-parameters.jl`: set fault parameters
    - `s04-solve.jl`: solve models
    - `scanfunc.jl`: functions for computing event catalogues from raw output data
    - `solve.job`: an example sbatch script for solving each model using one node (multi-threading)


- outputs (for catalogue only):

    All model outputs are in HDF5 file format. The names indicate: `catalogue-{parameter group}-{value}.h5`. LF denote left half fault, RF right half fault. Fields in the catalogue outputs are:

    - `t`: time step in second (downsampled)
    - `maxva`: max velocity in m/s (LF)
    - `maxvb`: max velocity in m/s (RF)
    - `mwa`: moment magnitude (LF)
    - `mwb`: moment magnitude (RF)
    - `ta`: event start time (LF)
    - `tb`: event start time (RF)
    - `sr`: seismic ratio
    - `ar`: afterslip ratio
    - `sync`: synchrony number


    All solution outputs (not included in this repo) contains three fields by default:

    - `t`: time step (s)
    - `v`: velocity (m/s)
    - `δ`: displacement (m)