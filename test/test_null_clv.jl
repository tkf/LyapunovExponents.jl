module TestNullCLV

using Base.Test
using Parameters: @with_kw, @unpack
using LyapunovExponents
using LyapunovExponents.CovariantVectors: is_finished
using LyapunovExponents: DEMOS, LEDemo, dimension

@with_kw struct NullCLVTest
    name::String
    prob::CLVProblem
    tolerance::Float64
    err_rate::Float64
end

NullCLVTest(demo::LEDemo;
            tolerance = nothing,
            err_rate = 0.0,
            kwargs...) = NullCLVTest(
    name = demo.example.name,
    prob = CLVProblem(demo.prob; kwargs...),
    tolerance = tolerance,
    err_rate = err_rate,
)

null_CLV_tests = [
    NullCLVTest(LyapunovExponents.lorenz_63();
                num_clv = 20000,
                num_forward_tran = 4000,
                num_backward_tran = 4000,
                err_rate = 0.002,  # 0.002 works with num_clv = 20000
                tolerance = 1e-2),
    # linz_sprott_99() requires a small `tspan` to make `tolerance`
    # smaller.
    NullCLVTest(LyapunovExponents.linz_sprott_99(tspan=0.1);
                num_clv = 20000,
                num_forward_tran = 4000,
                num_backward_tran = 4000,
                err_rate = 0.06,  # 0.06 works with num_clv = 2000000
                tolerance = 0.2),  # TODO: minimize
    NullCLVTest(LyapunovExponents.beer_95();
                num_clv = 200,
                num_forward_tran = 100,
                num_backward_tran = 100,
                err_rate = 0.01,  # 0.006 works with num_clv = 2000
                tolerance = 1e-2),
]

@time @testset "Null CLV: $(test.name)" for test in null_CLV_tests
    @unpack prob, tolerance, err_rate = test

    solver = init(prob; record=(:G, :C, :x))
    solve!(solver)

    dims = size(prob.Q0)
    x = solver.sol.x_history
    G = solver.sol.G_history
    C = solver.sol.C_history
    @test_broken length(x) == prob.num_clv
    @test_broken length(G) == prob.num_clv
    @test length(C) == prob.num_clv

    function ∂ₜ(u, prob = prob.phase_prob, t = 0.0)
        du = similar(u)
        prob.f(du, u, prob.p, t)
        return du
    end

    angles = collect(
        let ∂ₜxₙ = ∂ₜ(x[n]),
            vₙ = G[n] * C[n][:, 2]
            acos(abs(∂ₜxₙ' * vₙ) / norm(∂ₜxₙ))
        end
        for n in 1:prob.num_clv)

    @test mean(angles .> tolerance) <= err_rate

end

end
