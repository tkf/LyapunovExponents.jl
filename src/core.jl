import DifferentialEquations: init, solve, solve!, step!

"""
The Problem type represents the setup of the Lyapunov exponents calculation.
"""
abstract type AbstractLEProblem end

"""
The Relaxer type represents the calculation required for throwing away
the transient part of the dynamics.
"""
abstract type AbstractRelaxer end

"""
The Solver type represents the core Lyapunov Exponents (LE)
calculation.  The LE calculation is done by calling the in-place
mutating function [`solve!`](@ref).

The methods connecting the three principal types (Problem, Relaxer and
Solver) for the LE calculation are shown in the following diagram:

┌─ Problem ([`AbstractLEProblem`](@ref))                               \\\n
│     │                                                                \\\n
│     │ [`get_relaxer`](@ref), [`relaxed`](@ref)                       \\\n
│     ▼                                                                \\\n
│   Relaxer ([`AbstractRelaxer`](@ref)) ┄┄ ⟲ [`relax!`](@ref)        \\\n
│     │                                                                \\\n
│     │ [`LESolver`](@ref)                                             \\\n
│     │                                                                \\\n
│┄┄┄┄┄┄┄ [`init`](@ref), [`solve`](@ref)                         \\\n
│     │                                                                \\\n
│     ▼                                                                \\\n
└▶ Solver ([`AbstractLESolver`](@ref)) ┄┄ ⟲ [`solve!`](@ref),
                                                 [`step!`](@ref)

"""
abstract type AbstractLESolver end

dimension(prob) = length(prob.u0)

"""
    LEProblem(phase_prob;  <keyword arguments>)

# Arguments
- `phase_prob`: Phase space dynamics represented in the form of
  `ODEProblem` or `DiscreteProblem` from DifferentialEquations.jl.
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
struct LEProblem{DEP} <: AbstractLEProblem
    phase_prob::DEP
    num_tran
    dim_lyap
    Q0
    tangent_dynamics!

    function LEProblem(
            phase_prob::DEP;
            num_tran=1,
            dim_lyap=dimension(phase_prob),
            Q0=eye(dimension(phase_prob), dim_lyap),
            tangent_dynamics=nothing,
            ) where {DEP}
        new{DEP}(phase_prob, num_tran, dim_lyap, Q0, tangent_dynamics)
    end
end

struct Relaxer{LEP} <: AbstractRelaxer
    prob::LEP
    integrator
end

get_relaxer(prob::LEP; kwargs...) where {LEP <: AbstractLEProblem} =
    Relaxer{LEP}(prob, get_integrator(prob.phase_prob; kwargs...))

"""
    get_relaxer(prob::AbstractLEProblem; <keyword arguments>) :: AbstractRelaxer

Get a relaxer for a LE problem.
"""
function get_relaxer end

"""
    relax!(relaxer::AbstractRelaxer; <keyword arguments>)

Throwaway the transient part of the phase space dynamics of the LE
problem `prob`.
"""
function relax!(relaxer; progress=-1)
    integrator = relaxer.integrator
    num_tran = relaxer.prob.num_tran

    @showprogress_if(
        (progress >= 0), progress, "Transient dynamics...",
        for _ in 1:num_tran
            keepgoing!(integrator)
        end)

    relaxer
end

"""
    relaxed(prob::AbstractLEProblem; <keyword arguments>) :: AbstractRelaxer

Throwaway the transient part of the phase space dynamics of the LE
problem `prob`.

That is to say, convert a LE problem ([`AbstractLEProblem`](@ref)) to
a relaxer ([`AbstractRelaxer`](@ref)) and then call [`relax!`](@ref).
"""
relaxed(prob; progress=-1, kwargs...) =
    relax!(get_relaxer(prob; kwargs...); progress=progress)

function phase_tangent_state(relaxer::AbstractRelaxer)
    dim_lyap = relaxer.prob.dim_lyap
    Q0 = relaxer.prob.Q0
    x0 = last_state(relaxer)

    u0 = similar(x0, (length(x0), dim_lyap + 1))
    u0[:, 1] = x0
    u0[:, 2:end] = Q0
    u0
end

"""
    LESolver(integrator; <keyword arguments>)

A type representing the main calculation of Lyapunov Exponents (LE).
This struct holds all temporary state required for LE calculation.
"""
mutable struct LESolver{Intr} <: AbstractLESolver
    integrator::Intr
    exponents
    num_orth
    phase_state
    tangent_state
    sign_R

    function LESolver(
            integrator::Intr;
            phase_state=get_phase_state(integrator),
            tangent_state=get_tangent_state(integrator),
            dim_lyap=length(phase_state),
            ) where {Intr}
        num_orth = 0
        exponents = zeros(eltype(phase_state), dim_lyap)
        sign_R = Array{Bool}(dim_lyap)
        new{Intr}(
            integrator,
            exponents,
            num_orth,
            phase_state,
            tangent_state,
            sign_R,
        )
    end
end

function LESolver(prob::AbstractLEProblem, u0)
    phase_prob = prob.phase_prob

    de_prob_type = if isa(phase_prob, ODEProblem)
        ODEProblem
    elseif isa(phase_prob, DiscreteProblem)
        DiscreteProblem
    else
        error("Unknown problem type: $(typeof(phase_prob))")
    end
    # TODO: maybe use type parameter to get `de_prob_type`.

    tangent_dynamics! = prob.tangent_dynamics!
    if tangent_dynamics! == nothing
        phase_dynamics! = phase_prob.f
        tangent_dynamics! = PhaseTangentDynamics(phase_dynamics!, u0)
    end

    tangent_prob = de_prob_type(
        tangent_dynamics!,
        u0,
        phase_prob.tspan,
    )
    LESolver(
        tangent_prob,
    )
end

"""
    init(prob::AbstractLEProblem; <keyword arguments>) :: AbstractLESolver

Run phase space simulation to throw away the transient and then
construct a LE solver.
"""
init(prob::AbstractLEProblem; kwargs...) = LESolver(relaxed(prob; kwargs...))

@inline function keepgoing!(solver::AbstractLESolver)
    u0 = current_state(solver)
    u0[:, 1] = solver.phase_state
    u0[:, 2:end] = solver.tangent_state
    keepgoing!(solver.integrator, u0)
    u0 = current_state(solver)
    solver.phase_state[:] = @view u0[:, 1]
    solver.tangent_state[:] = @view u0[:, 2:end]
end

"""
S_n = ((n-1)/n) S_{n-1} + r_n / n
"""
@inline function lyap_add_R!(n, lyap, R)
    a = (n - 1) / n
    b = 1 - a
    for i = 1:length(lyap)
        lyap[i] = a * lyap[i] + b * log(R[i, i])
    end
end

""" A = A * diag(sgn) """
@inline function A_mul_sign!(A, sgn)
    for i = 1:size(A)[2]
        if !sgn[i]
            A[:, i] *= -1
        end
    end
    A
end

""" A = diag(sgn) * A """
@inline function sign_mul_A!(sgn, A)
    for i = 1:size(A)[1]
        if !sgn[i]
            A[i, :] *= -1
        end
    end
    A
end


""" sgn = sign(diag(A))  """
@inline function sign_diag!(sgn, A)
    for i = 1:size(A)[1]
        sgn[i] = A[i, i] >= 0
    end
    sgn
end

"""
    step!(solver::AbstractLESolver)

Evolve the dynamics and then do an orthonormalization.
"""
function step!(solver::AbstractLESolver)
    dim_lyap = length(solver.exponents)

    keepgoing!(solver)
    P = solver.tangent_state
    F = qrfact!(P)
    Q = F[:Q][:, 1:dim_lyap]
    R = F[:R]

    sign_R = solver.sign_R
    sign_diag!(sign_R, R)        # sign_R = diagm(sign(diag(R)))
    A_mul_sign!(Q, sign_R)       # Q = Q * sign_R
    sign_mul_A!(sign_R, R)       # R = signR * R

    n = (solver.num_orth += 1)
    lyap_add_R!(n, solver.exponents, R)

    solver.tangent_state = Q
end

"""
    solve!(solver::AbstractLESolver, num_attr; <keyword arguments>)

Do `num_attr` times of orthonormalization `step!(solver)`.
"""
function solve!(solver::AbstractLESolver, num_attr;
                progress=-1)
    @showprogress_if(
        (progress >= 0), progress, "Computing Lyapunov exponents...",
        for _ in 1:num_attr
            step!(solver)
        end)
    solver
end

"""
    solve(prob::AbstractLEProblem, num_attr; <keyword arguments>)
        :: AbstractLESolver

Initialize the solver ([`init`](@ref)) and then go through the LE
calculation ([`solve!`](@ref)).
"""
function solve(prob::AbstractLEProblem, num_attr; progress=-1, kwargs...)
    solver = init(prob; progress=progress, kwargs...)
    solve!(solver, num_attr; progress=progress)
end

"""
    lyapunov_exponents(solver)

Get the result of Lyapunov exponents calculation stored in `solver`.
"""
@inline function lyapunov_exponents(solver::AbstractLESolver)
    solver.exponents ./ t_chunk(solver)
end
# TODO: check if the dot here is meaningful (maybe define lyapunov_exponents!)

mutable struct LERecordingSolver <: AbstractLESolver
    solver
    exponents_history
    i_hist
end

function LERecordingSolver(solver, num_attr::Integer)
    dim_lyap = length(solver.exponents)
    exponents_history = similar(solver.exponents, (dim_lyap, num_attr))
    LERecordingSolver(solver, exponents_history, 0)
end

@inline function lyapunov_exponents(solver::LERecordingSolver)
    lyapunov_exponents(solver.solver)
end

function step!(solver::LERecordingSolver)
    step!(solver.solver)
    solver.i_hist += 1
    solver.exponents_history[:, solver.i_hist] .= lyapunov_exponents(solver)
end

function solve!(solver::LERecordingSolver; kwargs...)
    _, num_attr = size(solver.exponents_history)
    solver.i_hist = 0
    solve!(solver, num_attr; kwargs...)
end

Base.show(io::IO, solver::LERecordingSolver) = show(io, solver.solver)
Base.show(io::IO, solver::AbstractLESolver) =
    print(io,
          "#Orth.: ", solver.num_orth, ", ",
          "LEs: ", lyapunov_exponents(solver))
