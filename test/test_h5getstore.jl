using Test
using HDF5

@testset "property save and read" begin
    function test_equal(p)
        tmpfile = tempname()
        @store tmpfile p
        p′ = @getprop tmpfile
        @test p == p′
        rm(tmpfile)
    end

    ps = [
        SingleDofRSFProperty(rand(9)...),
        RateStateQuasiDynamicProperty([rand(9) for _ in 1: 4]..., rand(4)...),
        ]
    map(test_equal, ps)
end
