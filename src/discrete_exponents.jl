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

"""
    DiscreteLEProblem(phase_prob;  <keyword arguments>)

# Arguments
- `phase_prob::DiscreteProblem`: Phase space dynamics represented in the
  form of `DiscreteProblem` from DifferentialEquations.jl.
- `num_tran::Integer`: Number of iterations to through away to get rid
  of the transient dynamics.
- `dim_lyap::Integer`: Number of Lyapunov exponents to be calculated.
  Default to the full system dimension.
- `Q0::Array`: The initial guess of the Gram-Schmidt "Lyapunov vectors".
  Default to the identity matrix.
- `tangent_dynamics::Function`:
"""
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
