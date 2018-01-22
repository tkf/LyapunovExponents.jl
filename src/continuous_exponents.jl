using DifferentialEquations
using ForwardDiff

ODEIntegrator = OrdinaryDiffEq.ODEIntegrator

function get_integrator(prob; kwargs...)
    alg, extra_kwargs = default_algorithm(prob; kwargs...)
    init(prob, alg; kwargs..., extra_kwargs...)
end

"""
    ContinuousLEProblem(phase_prob; <keyword arguments>)

# Arguments
- `phase_prob::ODEProblem`: Phase space dynamics represented in the
  form of `ODEProblem` from DifferentialEquations.jl.
  `phase_prob.tspan` represents the inter-orthonormalization-interval.
- `num_tran::Integer`: Number of iterations to through away to get rid
  of the transient dynamics.
- `dim_lyap::Integer`: Number of Lyapunov exponents to be calculated.
  Default to the full system dimension.
- `Q0::Array`: The initial guess of the Gram-Schmidt "Lyapunov vectors".
  Default to the identity matrix.
- `tangent_dynamics::Function`: A vector field for solving phase space
   evolution *and* tangent space evolution together.  If this is not
   provided, `tangent_dynamics` is derived from `phase_prob.f`.  See
   also [`PhaseTangentDynamics`](@ref).
"""
const ContinuousLEProblem = LEProblem{<: ODEProblem}

"""
    ContinuousLEProblem(phase_prob; <keyword arguments>)
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
