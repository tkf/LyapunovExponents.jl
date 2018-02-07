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

    solver = init(prob)

    function ∂ₜ(u, prob = prob.phase_prob, t = 0.0)
        du = similar(u)
        prob.f(du, u, prob.p, t)
        return du
    end

    dims = size(prob.Q0)
    x = [Vector{Float64}(3) for _ in 1:prob.num_clv]
    G = [Matrix{Float64}(3, 3) for _ in 1:prob.num_clv]
    C = [Matrix{Float64}(3, 3) for _ in 1:prob.num_clv]

    forward = forward_dynamics!(solver)
    for (n, Gₙ) in zip(1:prob.num_clv, forward)
        x[n] .= phase_state(forward)
        G[n] .= Gₙ
    end
    @assert ! is_finished(forward)

    backward = backward_dynamics!(solver)
    @assert is_finished(forward)
    C[end] = CLV.C(backward)  # TODO: do this in the loop
    for (n, Cₙ) in indexed_backward_dynamics!(backward)
        C[n] .= Cₙ
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
