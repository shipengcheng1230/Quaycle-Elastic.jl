using Quaycle
using LinearAlgebra
using GmshTools
using Gmsh_SDK_jll

mf = gen_mesh(Val(:RectOkada), 80e3, 10e3, 80e3/512, 10e3/64, 90.0)
fmsh = joinpath(@__DIR__, splitext(basename(@__FILE__))[1] * "f.msh")
gen_gmsh_mesh(mf; filename=fmsh)