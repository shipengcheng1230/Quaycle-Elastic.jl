using Quaycle
using HDF5

gffile = joinpath(@__DIR__, "gf.h5")
include(joinpath(@__DIR__, "s01-domain.jl"))
@info "Fault elements: $(mf.nx * mf.nξ)"

ft = STRIKING()
λ = 3e10
μ = 3e10

@time ee = stress_greens_func(mf, λ, μ, ft; nrept=2, buffer_ratio=1)
isfile(gffile) && rm(gffile)
h5write(gffile, "ee", ee)