using DifferentialEquations
using ForwardDiff
using OrdinaryDiffEq: ODEIntegrator

function get_integrator(prob::ODEProblem; save_everystep=false, kwargs...)
    alg, extra_kwargs = default_algorithm(prob; kwargs...)
    return init(prob, alg; kwargs..., extra_kwargs...,
                save_everystep=save_everystep)
end
# See:
# - [[../../OrdinaryDiffEq/src/solve.jl::function init\b]]
# - http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html

const ContinuousLEProblem = LEProblem{ODEProblem}

"""
    ContinuousLEProblem(phase_dynamics!, u0 [, p];
                        t_attr=<number>, <keyword arguments>)

This is a short-hand notation for:

```julia
LEProblem(ODEProblem(phase_dynamics!, u0 [, p]), t_attr)
```

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
ContinuousLEProblem(phase_dynamics!, u0, p=nothing;
                    tspan=(0.0, 100.0), kwargs...) =
    ContinuousLEProblem(ODEProblem(phase_dynamics!, u0, tspan, p);
                        kwargs...)

init_phase_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 1]
init_tangent_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 2:end]

current_state(integrator::ODEIntegrator) = integrator.u
