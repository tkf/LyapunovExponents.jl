using DifferentialEquations
using ForwardDiff
using DiffEqBase: ODEProblem, SDEProblem, RODEProblem, DEIntegrator
using OrdinaryDiffEq: ODEIntegrator
using StochasticDiffEq: SDEIntegrator

""" Continuous-time problem. """
const CTProblem = Union{ODEProblem, SDEProblem, RODEProblem}
# TODO: maybe use Abstract*Problem; but note that
# AbstractDiscreteProblem <: AbstractODEProblem

function get_integrator(prob::CTProblem; save_everystep=false,
                        alg = nothing,
                        kwargs...)
    if alg === nothing
        alg, extra_kwargs = default_algorithm(prob; kwargs...)
    else
        extra_kwargs = []
    end
    return init(prob, alg; kwargs..., extra_kwargs...,
                save_everystep=save_everystep)
end
# See:
# - [[../../OrdinaryDiffEq/src/solve.jl::function init\b]]
# - http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html

const ContinuousLEProblem = LEProblem{ODEProblem}

"""
    ContinuousLEProblem(phase_dynamics, u0 [, p];
                        t_attr=<number>, <keyword arguments>)

This is a short-hand notation for:

```julia
LEProblem(ODEProblem(phase_dynamics, u0 [, p]), t_attr)
```

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
ContinuousLEProblem(phase_dynamics, u0, p=nothing;
                    tspan=(0.0, 100.0), kwargs...) =
    ContinuousLEProblem(ODEProblem(phase_dynamics, u0, tspan, p);
                        kwargs...)

init_phase_state(integrator::DEIntegrator) = integrator.sol.prob.u0[:, 1]
init_tangent_state(integrator::DEIntegrator) = integrator.sol.prob.u0[:, 2:end]

current_state(integrator::DEIntegrator) = integrator.u
