using Base.Test
using Base.Test: Fail, Error, Broken
using LyapunovExponents
using LyapunovExponents.TestTools: @test_isapprox_elemwise,
    @_test_isapprox_elemwise_result

@testset "@test_isapprox_elemwise" begin
    @test_isapprox_elemwise([0], [0])
    @test_isapprox_elemwise([1], [1.1], rtol=0.2)

    io = IOBuffer()
    result = @_test_isapprox_elemwise_result([0], [1], io=io)
    msg = readstring(seek(io, 0))
    @test contains(msg, "Not pairwise isapprox")
    @test result isa Fail

    # Marked as skip
    io = IOBuffer()
    result = @_test_isapprox_elemwise_result([0], [1], io=io, skip=true)
    @test result isa Broken

    # Marked as broken, but pass
    io = IOBuffer()
    result = @_test_isapprox_elemwise_result([0], [0], io=io, broken=true)
    @test result isa Error
    @test result.test_type == :test_unbroken

    # Marked as broken, and broken
    io = IOBuffer()
    result = @_test_isapprox_elemwise_result([0], [1], io=io, broken=true)
    @test result isa Broken
end
