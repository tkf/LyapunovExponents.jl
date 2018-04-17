using Base.Test
using LyapunovExponents
using LyapunovExponents: IntegratorError
using LyapunovExponents.TestTools: @test_nothrow

@testset "User interface" begin
    args = [
        (x...) -> nothing,  # f
        [0.0],              # u0
        nothing,            # p
    ]
    for lep in [DiscreteLEProblem, ContinuousLEProblem]
        @test_nothrow lep(args...; t_attr=1)
        try
            lep(args...)
            println("Calling $lep w/o t_attr didn't throw.")
            @test false
        catch err
            @test isa(err, ErrorException)
            @test contains(err.msg, "`t_attr` is required")
        end
    end
end

@testset "IntegratorError" begin
    prob = ContinuousLEProblem(
        (du, u, p, t) -> (@. du = exp(u)),  # f
        [0.0],    # u0 --- explosion time is 1.0
        nothing,  # p
        t_attr = 10.0,
        t_tran = 0,
    )
    err = try
        solve(prob, integrator_options=[:verbose => false])
        nothing
    catch err
        err
    end
    @test err isa IntegratorError
    @test err.retcode == :DtLessThanMin
end
