import DifferentialEquations: init, solve, solve!, step!

dimension(prob) = length(prob.u0)

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

function phase_tangent_state(relaxer::AbstractRelaxer)
    dim_lyap = relaxer.prob.dim_lyap
    Q0 = relaxer.prob.Q0
    x0 = last_state(relaxer)

    u0 = similar(x0, (length(x0), dim_lyap + 1))
    u0[:, 1] = x0
    u0[:, 2:end] = Q0
    u0
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
