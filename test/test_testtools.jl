using Base.Test
using LyapunovExponents
using LyapunovExponents.TestTools: @test_isapprox_elemwise

@testset "@test_isapprox_elemwise" begin
    @test_isapprox_elemwise([0], [0])
    @test_isapprox_elemwise([1], [1.1], rtol=0.2)
end
