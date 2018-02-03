# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
import DifferentialEquations: init, solve, solve!, step!

using DifferentialEquations: DEProblem

dimension(prob::DEProblem) = length(prob.u0)

"""
    get_relaxer(prob::AbstractLEProblem; <keyword arguments>) :: AbstractRelaxer

Get a relaxer for a LE problem.
"""
get_relaxer(prob::LEP; kwargs...) where {LEP <: AbstractLEProblem} =
    Relaxer{LEP}(prob, get_integrator(prob.phase_prob; kwargs...))

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

function phase_tangent_state(relaxer::Relaxer)
    x0 = last_state(relaxer)
    return phase_tangent_state(relaxer.prob, x0)
end

function phase_tangent_state(prob::LEProblem, x0 = prob.phase_prob.u0)
    dim_lyap = prob.dim_lyap
    Q0 = prob.Q0
    u0 = similar(x0, (length(x0), dim_lyap + 1))
    u0[:, 1] = x0
    u0[:, 2:end] = Q0
    u0
end

get_solver(relaxer::Relaxer; kwargs...) =
    get_solver(relaxer.prob, phase_tangent_state(relaxer); kwargs...)

function get_solver(prob::LEProblem{DEP},
                    u0 = phase_tangent_state(prob);
                    kwargs...) where {DEP}
    phase_prob = prob.phase_prob

    tangent_dynamics! = prob.tangent_dynamics!
    if tangent_dynamics! == nothing
        phase_dynamics! = phase_prob.f
        tangent_dynamics! = PhaseTangentDynamics(phase_dynamics!, u0)
    end

    tangent_prob = DEP(
        tangent_dynamics!,
        u0,
        phase_prob.tspan,
        phase_prob.p,
    )
    get_solver(tangent_prob, prob.num_attr; kwargs...)
end

function get_solver(tangent_prob::DEProblem, num_attr::Int;
                    record::Bool = false,
                    solver_type=nothing, kwargs...)
    if solver_type === nothing
        @assert ndims(tangent_prob.u0) == 2
        dim_lyap = size(tangent_prob.u0)[2] - 1
        if dim_lyap == 1
            solver_type = MLESolver
        else
            solver_type = LESolver
        end
    end

    solver = solver_type(get_integrator(tangent_prob), num_attr; kwargs...)
    if record
        solver = LERecordingSolver(solver)
    end
    return solver
end

"""
    init(prob::AbstractLEProblem; <keyword arguments>) :: AbstractLESolver
    init(relaxer::AbstractRelaxer; <keyword arguments>) :: AbstractLESolver

Run phase space simulation to throw away the transient and then
construct a LE solver.
"""
init(prob::AbstractLEProblem; progress = -1, kwargs...) =
    get_solver(relaxed(prob; progress = progress); kwargs...)
init(relaxer::AbstractRelaxer; kwargs...) = get_solver(relaxer; kwargs...)

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
    for i = 1:length(lyap)
        lyap[i] = lyap_add_r(n, lyap[i], R[i, i])
    end
end

@inline function lyap_add_r(n, lyap, r)
    return (n - 1) / n * lyap + log(r) / n
end

@inline function eye!(A)
    A .= 0
    for i in 1:min(size(A)...)
        @inbounds A[i, i] = 1
    end
    A
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
    keepgoing!(solver)
    post_evolve!(solver)
end

function post_evolve!(solver::LESolver)
    dim_lyap = length(solver.inst_exponents)
    P = solver.tangent_state
    F = qrfact!(P)
    Q = solver.Q
    A_mul_B!(F[:Q], eye!(Q))  # Q = Matrix(F[:Q])[...]; but lesser allocation
    R = solver.R = F[:R]

    sign_R = solver.sign_R
    sign_diag!(sign_R, R)        # sign_R = diagm(sign(diag(R)))
    A_mul_sign!(Q, sign_R)       # Q = Q * sign_R
    sign_mul_A!(sign_R, R)       # R = signR * R

    @assert size(R) == (dim_lyap, dim_lyap)
    for i in 1:dim_lyap
        @inbounds solver.inst_exponents[i] = log(R[i, i])
    end
    OnlineStats.fit!(solver.series, solver.inst_exponents)

    solver.num_orth += 1
    solver.tangent_state, solver.Q = Q, solver.tangent_state
end

function post_evolve!(solver::MLESolver)
    v = solver.tangent_state[:, 1]
    r = norm(v)
    n = (solver.num_orth += 1)
    solver.exponent = lyap_add_r(n, solver.exponent, r)
    solver.inst_exponent = log(r)  # TODO: don't re-calculate log(r)
    solver.tangent_state[:, 1] .= v ./ r
end

"""
    solve!(solver::AbstractLESolver; <keyword arguments>)

Do `solver.num_attr` times of orthonormalization `step!(solver)`.
"""
solve!(solver::AbstractLESolver; kwargs...) =
    forward!(solver, solver.num_attr; kwargs...)

function forward!(solver::AbstractLESolver, num_attr; progress=-1)
    @showprogress_if(
        (progress >= 0), progress, "Computing Lyapunov exponents...",
        for _ in 1:num_attr
            step!(solver)
        end)
    solver
end

"""
    solve(prob::AbstractLEProblem; <keyword arguments>)
        :: AbstractLESolver

Initialize the solver ([`init`](@ref)) and then go through the LE
calculation ([`solve!`](@ref)).
"""
function solve(prob::AbstractLEProblem;
               progress = -1,
               kwargs...)
    solver = init(prob; progress=progress, kwargs...)
    solve!(solver; progress=progress)
    return solver
end

"""
    lyapunov_exponents(solver)

Get the result of Lyapunov exponents calculation stored in `solver`.
"""
@inline function lyapunov_exponents(solver::LESolver)
    mean(solver.main_stat) ./ t_chunk(solver)
end

@inline function lyapunov_exponents(solver::MLESolver)
    [solver.exponent / t_chunk(solver)]
end

"""Get finite-time Lyapunov exponents (FTLE)"""
@inline ftle(solver::LESolver) = solver.inst_exponents ./ t_chunk(solver)
# TODO: check if the dot here is meaningful (maybe define ftle!)
# But probably this is not the performance-critical point.

@inline ftle(solver::MLESolver) = [solver.inst_exponent / t_chunk(solver)]
