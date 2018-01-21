using DifferentialEquations
using ForwardDiff

ODEIntegrator = OrdinaryDiffEq.ODEIntegrator

function get_integrator(prob; kwargs...)
    alg, extra_kwargs = default_algorithm(prob; kwargs...)
    init(prob, alg; kwargs..., extra_kwargs...)
end

struct ContinuousLEProblem <: AbstractLEProblem
    phase_prob
    num_tran
    dim_lyap
    Q0
    tangent_dynamics!

    function ContinuousLEProblem(
            phase_prob::ODEProblem;
            num_tran=1,
            dim_lyap=dimension(phase_prob),
            Q0=eye(dimension(phase_prob), dim_lyap),
            tangent_dynamics=nothing)
        new(phase_prob, num_tran, dim_lyap, Q0, tangent_dynamics)
    end
end

ContinuousLEProblem(phase_dynamics!, u0, tspan::Tuple; kwargs...) =
    ContinuousLEProblem(ODEProblem(phase_dynamics!, u0, tspan); kwargs...)

ContinuousLEProblem(phase_dynamics!, u0, tchunk::Real; kwargs...) =
    ContinuousLEProblem(phase_dynamics!, u0, (zero(tchunk), tchunk); kwargs...)

struct ContinuousRelaxer <: AbstractRelaxer
    prob
    integrator
end

get_relaxer(prob::ContinuousLEProblem; kwargs...) =
    ContinuousRelaxer(prob, get_integrator(prob.phase_prob; kwargs...))

last_state(relaxer::ContinuousRelaxer) = relaxer.integrator.sol[end]

get_phase_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 1]
get_tangent_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 2:end]
const ContinuousLESolver = LESolver{<: ODEIntegrator}

LESolver(tangent_prob::ODEProblem; kwargs...) =
    LESolver(get_integrator(tangent_prob); kwargs...)

LESolver(relaxer::ContinuousRelaxer) =
    LESolver(relaxer.prob, phase_tangent_state(relaxer))

"""
Continue solving the ODE problem from the last state.
"""
function keepgoing!(integrator::ODEIntegrator, u0=integrator.sol[end])
    reinit!(integrator, u0)
    solve!(integrator)
end
# TODO: Find out if this is a valid way for general algorithm

@inline function current_state(solver::ContinuousLESolver)
    solver.integrator.sol[end]
end

function t_chunk(solver::ContinuousLESolver)
    tspan = solver.integrator.sol.prob.tspan
    tspan[2] - tspan[1]
end
