using Test
using LinearAlgebra
using TensorOperations
using Quaycle:
    unit_dislocation,
    relative_velocity!,
    dτ_dt!,
    gen_alloc

@testset "Unit dislocation for plane fault types" begin
    @test unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end

@testset "2D FFT conv" begin
    mf = gen_mesh(Val(:RectOkada), 100.0, 100.0, 10.0, 10.0, 90.0)
    gf1 = stress_greens_func(mf, 3e10, 3e10, STRIKING(); fourier_domain=true)
    gf2 = stress_greens_func(mf, 3e10, 3e10, STRIKING(); fourier_domain=false)
    gf3 = Array{Float64}(undef, mf.nx, mf.nξ, mf.nx, mf.nξ)
    for l = 1: mf.nξ, k = 1: mf.nx, j = 1: mf.nξ, i = 1: mf.nx
        gf3[i,j,k,l] = gf2[abs(i-k)+1,j,l]
    end
    alloc = Quaycle.gen_alloc(gf1)
    v = rand(mf.nx, mf.nξ)
    vpl = 0.1
    relv = v .- vpl
    relative_velocity!(alloc, vpl, v)
    dτ_dt!(gf1, alloc)
    @tensor begin
        E[i,j] := gf3[i,j,k,l] * relv[k,l]
    end
    @test E ≈ alloc.dτ_dt
end
