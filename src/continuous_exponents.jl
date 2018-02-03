using DifferentialEquations
using ForwardDiff
using OrdinaryDiffEq: ODEIntegrator

function get_integrator(prob::ODEProblem; save_everystep=false, kwargs...)
    alg, extra_kwargs = default_algorithm(prob; save_everystep=save_everystep,
                                          kwargs...)
    init(prob, alg; kwargs..., extra_kwargs...)
end

const ContinuousLEProblem = LEProblem{ODEProblem}

"""
    ContinuousLEProblem(phase_dynamics!, u0, tspan [, p [, num_attr]];
                        <keyword arguments>)

Construct an `ODEProblem` and use it for `ContinuousLEProblem`.  If
`tspan` is a `Real` instead of a `Tuple`, then `(0, tspan)` is passed
as the `tspan` argument of `ODEProblem`.

For the list of usable keyword arguments, see [`LEProblem`](@ref).
"""
ContinuousLEProblem(phase_dynamics!, u0, tspan::Tuple, p=nothing,
                    args...; kwargs...) =
    ContinuousLEProblem(ODEProblem(phase_dynamics!, u0, tspan, p), args...;
                        kwargs...)

ContinuousLEProblem(phase_dynamics!, u0, tchunk::Real, args...; kwargs...) =
    ContinuousLEProblem(phase_dynamics!, u0, (zero(tchunk), tchunk), args...;
                        kwargs...)

const ContinuousRelaxer = Relaxer{<: ContinuousLEProblem}

last_state(relaxer::ContinuousRelaxer) = relaxer.integrator.sol[end]

init_phase_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 1]
init_tangent_state(integrator::ODEIntegrator) = integrator.sol.prob.u0[:, 2:end]
const ContinuousLESolver = AbstractLESolver{<: ODEIntegrator}

"""
Continue solving the ODE problem from the last state.
"""
function keepgoing!(integrator::ODEIntegrator, u0=integrator.sol[end])
    reinit!(integrator, u0)
    solve!(integrator)
end
# TODO: Find out if this is a valid way for general algorithm

current_state(integrator::ODEIntegrator) = integrator.sol[end]

@inline function current_state(solver::ContinuousLESolver)
    solver.integrator.sol[end]
end

function t_chunk(solver::ContinuousLESolver)
    tspan = solver.integrator.sol.prob.tspan
    tspan[2] - tspan[1]
end
