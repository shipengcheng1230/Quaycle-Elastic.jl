using Quaycle
using HDF5
using GmshTools
using Gmsh_SDK_jll

include(joinpath(@__DIR__, "s01-domain.jl"))

function set_mid_03(fname, x, type, val; avw=0.015, abvw=0.0045, Dc=7e-3, single=false, vwdepth=5e3, vwbarrier=false)
    cs = 3044.14
    vpl = 140e-3 / 365 / 86400 # 140 mm/yr
    v0 = 1e-6
    f0 = 0.6
    μ = 3e10
    λ = μ
    η = μ / 2cs
    ν = λ / 2(λ + μ)
    a = ones(mf.nx, mf.nξ) * avw # change
    b = ones(mf.nx, mf.nξ) * (avw - abvw) # change
    left_patch = (@. -x/2-20e3 ≤ mf.x ≤ -x/2) .* (@. (-1e3 - vwdepth) ≤ mf.z ≤ -1e3)'
    right_patch = (@. x/2 ≤ mf.x ≤ x/2+20e3) .* (@. (-1e3 - vwdepth) ≤ mf.z ≤ -1e3)'
    mid_patch = (@. -x/2 ≤ mf.x ≤ x/2) .* trues(mf.nξ)'
    bvw = avw + abvw
    b[left_patch] .= bvw # change
    if !single
        b[right_patch] .= bvw # change
    end
    if vwbarrier
        b[mid_patch] .= bvw
    end
    σmax = 50e6
    σ = [min(σmax, 1.5e6 + 18.0e3 * z) for z in -mf.z]
    σ = Matrix(repeat(σ, 1, mf.nx)')
    L = ones(mf.nx, mf.nξ) * Dc # change
    if type != nothing
        if type == "σ"
            for i in eachindex(σ)
                if mid_patch[i]
                    σ[i] = min(σ[i], val)
                end
            end
        elseif type == "a"
            a[mid_patch] .= val
        elseif type == "b"
            b[mid_patch] .= val
        elseif type == "L"
            L[mid_patch] .= val
        elseif type == "σc"
            σ[mid_patch] .= val
        end
    end
    pe = RateStateQuasiDynamicProperty(a, b, L, σ, η, vpl, f0, v0)
    f = joinpath(@__DIR__, fname)
    @store f pe
    cache = gmsh_vtk_output_cache(joinpath(@__DIR__, "dom01f.msh"), mf)
    vtk_output(joinpath(@__DIR__, splitext(f)[1]), pe, cache)
end

begin # for smaller VW
    set_mid_03("pf-s-w10.h5", 10e3, nothing, 1.0, single=true,   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-00.h5",    10e3, nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s200.h5",  10e3, "σc", 200e6,                 abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s150.h5",  10e3, "σc", 150e6,                 abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s100.h5",  10e3, "σc", 100e6,                 abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s80.h5",   10e3, "σc", 80e6,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s60.h5",   10e3, "σc", 60e6,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s25.h5",   10e3, "σ", 25e6,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s10.h5",   10e3, "σ", 10e6,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s5.h5",    10e3, "σ", 5e6,                    abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-s1.h5",    10e3, "σ", 1e6,                    abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w2.h5",    2e3,  nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w3.h5",    3e3,  nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w5.h5",    5e3,  nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w15.h5",   15e3, nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w20.h5",   20e3, nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w30.h5",   30e3, nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-w40.h5",   40e3, nothing, 1.0,                abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-b4.0.h5",  10e3, "b", 0.015 - 0.0035 * 4.0,   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-b3.0.h5",  10e3, "b", 0.015 - 0.0035 * 3.0,   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-b2.0.h5",  10e3, "b", 0.015 - 0.0035 * 2.0,   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-b1.35.h5", 10e3, "b", 0.015 - 0.0035 * 1.35,  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-b0.5.h5",  10e3, "b", 0.015 - 0.0035 * 0.5,   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-b0.2.h5",  10e3, "b", 0.015 - 0.0035 * 0.2,   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-a4.0.h5",  10e3, "a", 0.0115 + 0.0035 * 4.0,  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-a3.0.h5",  10e3, "a", 0.0115 + 0.0035 * 3.0,  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-a2.0.h5",  10e3, "a", 0.0115 + 0.0035 * 2.0,  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-a1.35.h5", 10e3, "a", 0.0115 + 0.0035 * 1.35, abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-a0.5.h5",  10e3, "a", 0.0115 + 0.0035 * 0.5,  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-a0.2.h5",  10e3, "a", 0.0115 + 0.0035 * 0.2,  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l1.0.h5",  10e3, "L", 1e-3,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l2.0.h5",  10e3, "L", 2e-3,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l5.0.h5",  10e3, "L", 5e-3,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l7.0.h5",  10e3, "L", 7e-3,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l8.0.h5",  10e3, "L", 8e-3,                   abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l12.0.h5", 10e3, "L", 12e-3,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l16.0.h5", 10e3, "L", 16e-3,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l20.0.h5", 10e3, "L", 20e-3,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l40.0.h5", 10e3, "L", 40e-3,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l60.0.h5", 10e3, "L", 60e-3,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
    set_mid_03("pf-l80.0.h5", 10e3, "L", 80e-3,                  abvw=0.0035, Dc=4e-3, avw=0.015, vwdepth=5e3)
end