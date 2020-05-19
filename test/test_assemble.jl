using Test

@testset "Prob Assemble" begin
    @testset "1D fault" begin
        mesh = gen_mesh(Val(:LineOkada), 10., 2.0, 45.0)
        gf = stress_greens_func(mesh, 1.0, 1.0, DIPPING())
        p = RateStateQuasiDynamicProperty([rand(mesh.nξ) for _ in 1: 4]..., rand(4)...)
        u0 = ArrayPartition([rand(mesh.nξ) for _ in 1: 3]...)
        prob = assemble(gf, p, u0, (0., 1.0); flf=CForm())
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

    @testset "2D fault" begin
        mesh = gen_mesh(Val(:RectOkada), 10., 10., 2., 2., 90.)
        gf = stress_greens_func(mesh, 1.0, 1.0, STRIKING(); buffer_ratio=1.0)
        p = RateStateQuasiDynamicProperty([rand(mesh.nx, mesh.nξ) for _ in 1: 4]..., rand(4)...)
        u0 = ArrayPartition([rand(mesh.nx, mesh.nξ) for _ in 1: 3]...)
        prob = assemble(gf, p, u0, (0., 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end
end
