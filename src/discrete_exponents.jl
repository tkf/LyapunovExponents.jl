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

function keepgoing!(diter::DiscreteIterator, u0=diter.u0)
    tmin, tmax = diter.prob.tspan
    f = diter.prob.f
    u1 = diter.u1
    for t in tmin:tmax
        f(t, u0, u1)
        u0, u1 = u1, u0
    end
    diter.u0 = u0
    diter.u1 = u1
end

const DiscreteLEProblem = LEProblem{<: DiscreteProblem}

"""
    DiscreteLEProblem(phase_dynamics!, u0, tspan; <keyword arguments>)

Construct a `DiscretProblem` and use it for `DiscreteLEProblem`.

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
DiscreteLEProblem(phase_dynamics!, u0, tspan; kwargs...) =
    LEProblem(DiscreteProblem(phase_dynamics!, u0, tspan); kwargs...)

struct DiscreteRelaxer <: AbstractRelaxer
    prob
    integrator
end

get_relaxer(prob::DiscreteLEProblem) =
    DiscreteRelaxer(prob, DiscreteIterator(prob.phase_prob))

last_state(relaxer::DiscreteRelaxer) = relaxer.integrator.u0

get_phase_state(integrator::DiscreteIterator) = integrator.u0[:, 1]
get_tangent_state(integrator::DiscreteIterator) = integrator.u0[:, 2:end]
const DiscreteLESolver = LESolver{<: DiscreteIterator}

LESolver(tangent_prob::DiscreteProblem; kwargs...) =
    LESolver(DiscreteIterator(tangent_prob); kwargs...)

LESolver(relaxer::DiscreteRelaxer; kwargs...) =
    LESolver(relaxer.prob, phase_tangent_state(relaxer))

@inline function current_state(solver::DiscreteLESolver)
    solver.integrator.u0
end

function t_chunk(solver::DiscreteLESolver)
    tspan = solver.integrator.prob.tspan
    tspan[2] - tspan[1] + 1
end
