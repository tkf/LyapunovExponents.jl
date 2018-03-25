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

DiffEqBase.set_u!(diter::DiscreteIterator, u0) = diter.u0 = u0

function step!(diter::DiscreteIterator, dt::Int = 1, _ignored::Bool = false)
    f = diter.prob.f
    p = diter.prob.p
    u0 = diter.u0
    u1 = diter.u1
    for t in 0:dt-1
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
    DiscreteLEProblem(phase_dynamics, u0 [, p];
                      t_attr=<number>, <keyword arguments>)

This is a short-hand notation for:

```julia
LEProblem(DiscreteProblem(phase_dynamics, u0 [, p]), t_attr)
```

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
DiscreteLEProblem(phase_dynamics, u0, p=nothing;
                  tspan=(0, typemax(Int)), kwargs...) =
    DiscreteLEProblem(DiscreteProblem(phase_dynamics, u0, tspan, p);
                      kwargs...)

init_phase_state(integrator::DiscreteIterator) = integrator.u0[:, 1]
init_tangent_state(integrator::DiscreteIterator) = integrator.u0[:, 2:end]

current_state(integrator::DiscreteIterator) = integrator.u0
