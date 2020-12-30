using Quaycle
using HDF5
using LinearAlgebra

# mkl performs better than openblas if on Intel CPUs
if BLAS.vendor() == :openblas64
    @info ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
end

function solve_from_para_01(pf, output; basedir="", offsetinit=2.5, stride=100, yearto=600.0)
    @info "Loading $(pf) ..."
    # I didn't check if the folder exists, change it to your like
    pff = joinpath(dirname(@__DIR__), "para", pf)
    pe = @getprop pff

    gffile = joinpath(@__DIR__, "gf.h5")
    @info "Loading Green's function " * basename(gffile) * " ..."
    gf = h5read(gffile, "ee")

    @info "Initializing ..."
    if isa(offsetinit, Number)
        @info "Manual offset $(offsetinit) ..."
        vinit = pe.vpl .* ones(size(pe.a))
        θinit = pe.L ./ vinit
        halflen = size(pe.a, 1) ÷ 2
        θinit[1:halflen,:] ./= 1.1
        θinit[halflen+1:end,:] ./= offsetinit # for smaller VW
    elseif isa(offsetinit, String)
        @info "Offset from icfile $(offsetinit) ..."
        icfile = joinpath(dirname(@__DIR__), "out", offsetinit)
        vinit = h5read(icfile, "v")
        θinit = h5read(icfile, "θ")
    end
    δinit = zeros(size(pe.a))
    uinit = ArrayPartition(vinit, θinit, δinit)
    prob = assemble(gf, pe, uinit, (0.0, yearto * 365 * 86400))
    @info "Solving $(output) ..."
    # I didn't check if the folder exists, change it to your like
    output_ = joinpath(dirname(@__DIR__), "out", basedir, output)
    @time sol = wsolve(prob, VCABM5(), output_, 1000, 𝐕𝚯𝚫, ["v", "θ", "δ"], "t";
        reltol=1e-6, abstol=1e-8, dtmax=0.2*365*86400, dt=1e-8, maxiters=1e9, stride=stride, force=false)
end

include(joinpath(@__DIR__, "scanfunc.jl"))

model2output_orgin = Dict(
    "otf-s-s25.h5" =>   ( "pf-s-s25.h5",),
    "otf-s-s60.h5" =>   ( "pf-s-s60.h5", ),
    "otf-s-a0.5.h5" =>  ( "pf-s-a0.5.h5", ),

    "otf-00.h5"     => ( "pf-00.h5",  ),

    "otf-s200.h5"   => ( "pf-s200.h5",),
    "otf-s150.h5"   => ( "pf-s150.h5",),
    "otf-s100.h5"   => ( "pf-s100.h5",),
    "otf-s80.h5"    => ( "pf-s80.h5", ),
    "otf-s60.h5"    => ( "pf-s60.h5", ),
    "otf-s25.h5"    => ( "pf-s25.h5", ),
    "otf-s10.h5"    => ( "pf-s10.h5", ),
    "otf-s5.h5"     => ( "pf-s5.h5",  ),
    "otf-s1.h5"     => ( "pf-s1.h5",  ),

    "otf-w2.h5"     => ( "pf-w2.h5",  ),
    "otf-w3.h5"     => ( "pf-w3.h5",  ),
    "otf-w5.h5"     => ( "pf-w5.h5",  ),
    "otf-w15.h5"    => ( "pf-w15.h5", ),
    "otf-w20.h5"    => ( "pf-w20.h5", ),
    "otf-w30.h5"    => ( "pf-w30.h5", ),
    "otf-w40.h5"    => ( "pf-w40.h5", ),

    "otf-b4.0.h5"   => ( "pf-b4.0.h5",),
    "otf-b3.0.h5"   => ( "pf-b3.0.h5",),
    "otf-b2.0.h5"   => ( "pf-b2.0.h5",),
    "otf-b1.35.h5"  => ( "pf-b1.35.h5"),
    "otf-b0.5.h5"   => ( "pf-b0.5.h5",),
    "otf-b0.2.h5"   => ( "pf-b0.2.h5",),

    "otf-a4.0.h5"   => ( "pf-a4.0.h5",),
    "otf-a3.0.h5"   => ( "pf-a3.0.h5",),
    "otf-a2.0.h5"   => ( "pf-a2.0.h5",),
    "otf-a1.35.h5"  => ( "pf-a1.35.h5"),
    "otf-a0.5.h5"   => ( "pf-a0.5.h5",),
    "otf-a0.2.h5"   => ( "pf-a0.2.h5",),

    "otf-l1.0.h5"   => ( "pf-l1.0.h5",),
    "otf-l2.0.h5"   => ( "pf-l2.0.h5",),
    "otf-l5.0.h5"   => ( "pf-l5.0.h5",),
    "otf-l7.0.h5"   => ( "pf-l7.0.h5",),
    "otf-l8.0.h5"   => ( "pf-l8.0.h5",),
    "otf-l12.0.h5"  => ( "pf-l12.0.h5"),
    "otf-l16.0.h5"  => ( "pf-l16.0.h5"),
    "otf-l20.0.h5"  => ( "pf-l20.0.h5"),
    "otf-l40.0.h5"  => ( "pf-l40.0.h5"),
    "otf-l60.0.h5"  => ( "pf-l60.0.h5"),
    "otf-l80.0.h5"  => ( "pf-l80.0.h5"),
)

mkey = keys(model2output_orgin) |> collect |> sort
k = mkey[parse(Int64, ARGS[1])] # output name
v = model2output_orgin[k] # prop name

# basedir = "ic2.5" | "ic1.5" | "ic1.1"
# offsetinit = 2.5 | 1.5 | 1.1
basedir = "ic2.5"
offsetinit = 2.5
solve_from_para_01(v[1], k, yearto=600.0, stride=100, offsetinit=offsetinit, basedir=basedir)

# edge = 5e3 specifically for VS_sigma = 1 MPa
scan_output_01(k; basedir=basedir, edge=1e3)
slip_ratio(k, -6e3)
