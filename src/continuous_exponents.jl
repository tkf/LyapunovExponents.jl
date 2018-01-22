using DifferentialEquations
using ForwardDiff

ODEIntegrator = OrdinaryDiffEq.ODEIntegrator

function get_integrator(prob; kwargs...)
    alg, extra_kwargs = default_algorithm(prob; kwargs...)
    init(prob, alg; kwargs..., extra_kwargs...)
end

const ContinuousLEProblem = LEProblem{<: ODEProblem}

"""
    ContinuousLEProblem(phase_dynamics!, u0, tspan; <keyword arguments>)

Construct an `ODEProblem` and use it for `ContinuousLEProblem`.  If
`tspan` is a `Real` instead of a `Tuple`, then `(0, tspan)` is passed
as the `tspan` argument of `ODEProblem`.

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
ContinuousLEProblem(phase_dynamics!, u0, tspan::Tuple; kwargs...) =
    LEProblem(ODEProblem(phase_dynamics!, u0, tspan); kwargs...)

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
