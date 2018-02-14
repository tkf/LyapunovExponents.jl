using DifferentialEquations
using ForwardDiff
using OrdinaryDiffEq: ODEIntegrator

function get_integrator(prob::ODEProblem; save_everystep=false, kwargs...)
    alg, extra_kwargs = default_algorithm(prob; save_everystep=save_everystep,
                                          kwargs...)
    init(prob, alg; kwargs..., extra_kwargs...)
end
# See:
# - [[../../OrdinaryDiffEq/src/solve.jl::function init\b]]
# - http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html

const ContinuousLEProblem = LEProblem{ODEProblem}

"""
    ContinuousLEProblem(phase_dynamics!, u0, tspan [, p [, num_attr]];
                        <keyword arguments>)

This is a short-hand notation for:

```julia
LEProblem(ODEProblem(...) [, num_attr]; ...)
```

If `tspan` is a `Real` instead of a `Tuple`, then `(0, tspan)` is
passed as the `tspan` argument of `ODEProblem`.

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
ContinuousLEProblem(phase_dynamics!, u0, tspan::Tuple, p=nothing,
                    args...; kwargs...) =
    ContinuousLEProblem(ODEProblem(phase_dynamics!, u0, tspan, p), args...;
                        kwargs...)

ContinuousLEProblem(phase_dynamics!, u0, tchunk::Real, args...; kwargs...) =
    ContinuousLEProblem(phase_dynamics!, u0, (zero(tchunk), tchunk), args...;
                        kwargs...)

init_phase_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 1]
init_tangent_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 2:end]

"""
Continue solving the ODE problem from the last state.
"""
function keepgoing!(integrator::ODEIntegrator, u0=integrator.sol[end])
    reinit!(integrator, u0)
    solve!(integrator)
end
# TODO: Find out if this is a valid way for general algorithm

current_state(integrator::ODEIntegrator) = integrator.sol[end]
get_tspan(integrator::ODEIntegrator) = integrator.sol.prob.tspan
