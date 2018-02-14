mutable struct DiscreteIterator
    prob
    u0
    u1

    function DiscreteIterator(prob)
        u0 = copy(prob.u0)
        u1 = similar(prob.u0)
        new(prob, u0, u1)
    end
end

get_integrator(prob::DiscreteProblem) = DiscreteIterator(prob)

function keepgoing!(diter::DiscreteIterator, u0=diter.u0)
    tmin, tmax = diter.prob.tspan
    f = diter.prob.f
    p = diter.prob.p
    u1 = diter.u1
    for t in tmin:tmax-1
        f(u1, u0, p, t)
        u0, u1 = u1, u0
    end
    diter.u0 = u0
    diter.u1 = u1
end
# See:
# - [[../../OrdinaryDiffEq/src/perform_step/fixed_timestep_perform_step.jl::perform_step!.*DiscreteConstantCache]]
# - [[../../OrdinaryDiffEq/src/solve.jl::solve!.*::ODEIntegrator]]
# - http://devdocs.juliadiffeq.org/latest/contributing/diffeq_internals.html

const DiscreteLEProblem = LEProblem{DiscreteProblem}

"""
    DiscreteLEProblem(phase_dynamics!, u0, tspan [, p [, num_attr]];
                      <keyword arguments>)

This is a short-hand notation for:

```julia
LEProblem(DiscreteProblem(...) [, num_attr]; ...)
```

If `tspan` is a `Integer` instead of a `Tuple`, then `(0, tspan)` is
passed as the `tspan` argument of `DiscreteProblem`.

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
DiscreteLEProblem(phase_dynamics!, u0, tspan::Tuple, p=nothing,
                  args...; kwargs...) =
    DiscreteLEProblem(DiscreteProblem(phase_dynamics!, u0, tspan, p),
                      args...; kwargs...)

DiscreteLEProblem(phase_dynamics!, u0, tchunk::Integer, args...; kwargs...) =
    DiscreteLEProblem(phase_dynamics!, u0, (zero(tchunk), tchunk), args...;
                      kwargs...)

init_phase_state(integrator::DiscreteIterator) = integrator.u0[:, 1]
init_tangent_state(integrator::DiscreteIterator) = integrator.u0[:, 2:end]

current_state(integrator::DiscreteIterator) = integrator.u0
get_tspan(integrator::DiscreteIterator) = integrator.prob.tspan
