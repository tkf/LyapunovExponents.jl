"""
Test that `DiscreteIterator` has the same semantics as
`DifferentialEquations.solve`.
"""
module TestDiscreteIterator

using Base.Test
using OrdinaryDiffEq: DiscreteProblem
using LyapunovExponents: DiscreteIterator, keepgoing!, change_tspan,
    solve_discrete

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

    sol = solve_discrete(prob)

    @test iter.u0 ≈ sol[end]
end

end
