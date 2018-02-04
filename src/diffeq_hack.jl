using OrdinaryDiffEq: DiscreteProblem, FunctionMap, solve

change_tspan(prob, tspan::Tuple) =
    DiscreteProblem(prob.f, prob.u0, tspan, prob.p)

change_tspan(prob, tchunk::Integer) = change_tspan(prob, (0, tchunk))

"""Convert `prob.tspan` to `Tuple{Float64, Float64}`."""
float_tspan(prob) = change_tspan(prob, Tuple{Float64, Float64}(prob.tspan))

"""
    solve_discrete(prob::DiscreteProblem; <keyword arguments>)

Solve prob using DifferentialEquations.jl with a bit of hack to make
it work.
See: https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/251
"""
solve_discrete(prob; kwargs...) =
    solve(
        # Convert tspan to Tuple{Float64, Float64}:
        float_tspan(prob),
        # Avoid DifferentialEquations.jl to repeat the initial condition:
        FunctionMap(scale_by_time=false);
        kwargs...)
