"""
Test that `DiscreteIterator` has the same semantics as
`DifferentialEquations.solve`.
"""
module TestDiscreteIterator

using Base.Test
using OrdinaryDiffEq: DiscreteProblem, FunctionMap, solve
using LyapunovExponents: DiscreteIterator, keepgoing!

change_tspan(prob, tspan::Tuple) =
    DiscreteProblem(prob.f, prob.u0, tspan, prob.p)

change_tspan(prob, tchunk::Integer) = change_tspan(prob, (0, tchunk))

one_step() = DiscreteProblem(
    (u_next, u, p, t) -> u_next .= 2 .* u,
    [1.0],
    (0, 1),
)

two_steps() = change_tspan(one_step(), 2)
ten_steps() = change_tspan(one_step(), 10)


@testset "$make_prob" for make_prob in [
            one_step,
            two_steps,
            ten_steps,
        ]
    prob = make_prob()

    iter = DiscreteIterator(prob)
    keepgoing!(iter)

    # Solve prob using DifferentialEquations.jl.
    # A bit of hack is required to make it work.
    # See: https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/251
    sol = solve(
        # Convert tspan to Tuple{Float64, Float64}:
        change_tspan(prob, Tuple{Float64, Float64}(prob.tspan)),
        # Avoid DifferentialEquations.jl to repeat the initial condition:
        FunctionMap(scale_by_time=false),
    )

    @test iter.u0 ≈ sol[end]
end

end
