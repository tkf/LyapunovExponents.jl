using Base.Test
using Base.Test: Fail, Error, Broken
using LyapunovExponents
using LyapunovExponents: objname, dimension
using LyapunovExponents.TestTools: test_same_dynamics,
    test_tangent_dynamics_against_autodiff,
    @test_isapprox_elemwise, @_test_isapprox_elemwise_result

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

@testset "test_same_dynamics: $(objname(f))" for f in [
        LyapunovExponents.henon_map,
        LyapunovExponents.van_der_pol,
        ]
    p1 = f().prob.phase_prob
    p2 = f().prob.phase_prob
    num_u = 5  # number of random states at which the systems are tested
    test_same_dynamics(p1, p2, num_u)
    test_same_dynamics(p1, p2, num_u; evolve=true)
    test_same_dynamics(p1, p2, [p1.u0, p1.u0 + 0.1])
end

@testset "test_tangent_dynamics_against_autodiff: dim_lyap=$dim_lyap" for
        dim_lyap in 1:3
    prob = LEProblem(LyapunovExponents.lorenz_63().example; dim_lyap=dim_lyap)
    @test size(prob.tangent_prob.u0) == (3, 1 + dim_lyap)
    num_u = 5
    test_tangent_dynamics_against_autodiff(prob, num_u)
end
