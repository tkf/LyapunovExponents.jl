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

struct DiscreteLEProblem <: AbstractLEProblem
    phase_prob
    num_tran
    dim_lyap
    Q0
    tangent_dynamics!

    function DiscreteLEProblem(
            phase_prob::DiscreteProblem;
            num_tran=1,
            dim_lyap=dimension(phase_prob),
            Q0=eye(dimension(phase_prob), dim_lyap),
            tangent_dynamics=nothing
            )
        new(phase_prob, num_tran, dim_lyap, Q0, tangent_dynamics)
    end
end

DiscreteLEProblem(phase_dynamics!, u0, tspan; kwargs...) =
    DiscreteLEProblem(DiscreteProblem(phase_dynamics!, u0, tspan); kwargs...)

struct DiscreteRelaxer <: AbstractRelaxer
    prob
    integrator
end

get_relaxer(prob::DiscreteLEProblem) =
    DiscreteRelaxer(prob, DiscreteIterator(prob.phase_prob))

last_state(relaxer::DiscreteRelaxer) = relaxer.integrator.u0

mutable struct DiscreteLESolver <: AbstractLESolver
    integrator
    exponents
    num_orth
    phase_state
    tangent_state
    sign_R

    function DiscreteLESolver(
            integrator::DiscreteIterator;
            phase_state=integrator.u0[:, 1],
            tangent_state=integrator.u0[:, 2:end],
            dim_lyap=length(phase_state))
        num_orth = 0
        exponents = zeros(eltype(phase_state), dim_lyap)
        sign_R = Array{Bool}(dim_lyap)
        new(
            integrator,
            exponents,
            num_orth,
            phase_state,
            tangent_state,
            sign_R,
        )
    end
end

DiscreteLESolver(tangent_prob::DiscreteProblem; kwargs...) =
    DiscreteLESolver(DiscreteIterator(tangent_prob); kwargs...)

function DiscreteLESolver(prob::DiscreteLEProblem, u0)
    phase_prob = prob.phase_prob

    tangent_dynamics! = prob.tangent_dynamics!
    if tangent_dynamics! == nothing
        phase_dynamics! = phase_prob.f
        tangent_dynamics! = PhaseTangentDynamics(phase_dynamics!, u0)
    end

    tangent_prob = DiscreteProblem(
        tangent_dynamics!,
        u0,
        phase_prob.tspan,
    )
    DiscreteLESolver(
        tangent_prob,
    )
end

function DiscreteLESolver(relaxer::DiscreteRelaxer; kwargs...)
    DiscreteLESolver(relaxer.prob, phase_tangent_state(relaxer))
end

function init(prob::DiscreteLEProblem; kwargs...)
    DiscreteLESolver(relaxed(prob; kwargs...))
end

@inline function current_state(solver::DiscreteLESolver)
    solver.integrator.u0
end

function t_chunk(solver::DiscreteLESolver)
    tspan = solver.integrator.prob.tspan
    tspan[2] - tspan[1] + 1
end
