using DifferentialEquations: DEProblem
using Distributions: TDist, Normal, cquantile

dimension(prob::DEProblem) = length(prob.u0)

PhaseRelaxer(prob::LEProblem,
             ::LEProblem,
             ::LESolution) =
    PhaseRelaxer(get_integrator(prob.phase_prob),
                 prob.num_tran)

step!(relaxer::PhaseRelaxer) = keepgoing!(relaxer.integrator)
Base.length(relaxer::PhaseRelaxer) = relaxer.num_tran

function phase_tangent_state(prob::LEProblem, x0 = prob.phase_prob.u0)
    dim_lyap = prob.dim_lyap
    Q0 = prob.Q0
    u0 = similar(x0, (length(x0), dim_lyap + 1))
    u0[:, 1] = x0
    u0[:, 2:end] = Q0
    u0
end

function get_tangent_dynamics(prob, u0 = phase_tangent_state(prob))
    phase_prob = prob.phase_prob
    tangent_dynamics! = prob.tangent_dynamics!
    if tangent_dynamics! == nothing
        phase_dynamics! = phase_prob.f
        tangent_dynamics! = PhaseTangentDynamics(phase_dynamics!, u0)
    end
    return tangent_dynamics!
end

function get_tangent_prob(prob::LEProblem{DEP},
                          u0 = phase_tangent_state(prob)) where {DEP}
    phase_prob = prob.phase_prob
    return DEP(
        get_tangent_dynamics(prob),
        u0,
        phase_prob.tspan,
        phase_prob.p,
    )
end

current_state(relaxer::Union{PhaseRelaxer, AbstractRenormalizer}) =
    current_state(relaxer.integrator)

function get_tangent_integrator(prob::LEProblem, relaxer)
    u0 = phase_tangent_state(prob, current_state(relaxer))
    tangent_prob = get_tangent_prob(prob, u0)
    return get_integrator(tangent_prob)
end

TangentRenormalizer(relaxer::PhaseRelaxer, prob::LEProblem, sol::LESolution) =
    TangentRenormalizer(get_tangent_integrator(prob, relaxer),
                        prob.num_attr, sol)

MLERenormalizer(relaxer::PhaseRelaxer, prob::LEProblem, sol::LESolution) =
    MLERenormalizer(get_tangent_integrator(prob, relaxer), sol, prob.num_attr)

Base.length(stage::AbstractRenormalizer) = stage.num_attr

@inline function keepgoing!(solver::AbstractRenormalizer)
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
    step!(solver::AbstractRenormalizer)

Evolve the dynamics and then do an orthonormalization.
"""
function step!(solver::AbstractRenormalizer)
    keepgoing!(solver)
    post_evolve!(solver)
    stage = solver
    record!(stage, Val{:OS})
    record!(stage, Val{:FTLE})
end

function post_evolve!(solver::TangentRenormalizer)
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

    solver.i += 1
    solver.tangent_state, solver.Q = Q, solver.tangent_state
    # At this point:
    # solver.tangent_state === Q == 𝑮ₙ₊ₖ
    # solver.Q is a mess, as qrfact! use it as a "buffer"
end

function post_evolve!(solver::MLERenormalizer)
    v = solver.tangent_state[:, 1]
    r = norm(v)
    n = (solver.i += 1)
    solver.exponent = lyap_add_r(n, solver.exponent, r)
    stage = solver
    stage.inst_exponent = log(r) / t_chunk(stage)
    # TODO: don't re-calculate log(r)
    solver.tangent_state[:, 1] .= v ./ r
end

function t_chunk(stage)
    tspan = get_tspan(stage.integrator)
    return tspan[2] - tspan[1]
end

function record!(stage::TangentRenormalizer{<: LESolRecOS}, ::Type{Val{:OS}})
    dim_lyap = length(stage.inst_exponents)
    R = stage.R
    dt = t_chunk(stage)
    @assert size(R) == (dim_lyap, dim_lyap)
    for i in 1:dim_lyap
        @inbounds stage.inst_exponents[i] = log(R[i, i]) / dt
    end
    OnlineStats.fit!(stage.sol.series, stage.inst_exponents)
end

function record!(stage::MLERenormalizer{<: LESolRecOS}, ::Type{Val{:OS}})
    # FIXME: Assuming stage.sol.main is VecMean and alike
    OnlineStats.fit!(stage.sol.series, [stage.inst_exponent])
end

record!(stage::AbstractRenormalizer{<: LESolRecFTLE}, ::Type{Val{:FTLE}}) =
    stage.sol.ftle_history[stage.i] .= ftle(stage)

"""
    lyapunov_exponents(solver)

Get the result of Lyapunov exponents calculation stored in `solver`.
"""
lyapunov_exponents(solver::Union{LESolverRecOS,
                                 AbstractRenormalizer{<: LESolRecOS}}) =
    lyapunov_exponents(solver.sol)

lyapunov_exponents(sol::LESolRecOS) = mean(sol.main_stat)

@inline function lyapunov_exponents(solver::MLERenormalizer)
    [solver.exponent / t_chunk(solver)]
end

"""Get finite-time Lyapunov exponents (FTLE)"""
ftle(stage::TangentRenormalizer) = stage.inst_exponents
ftle(stage::MLERenormalizer) = [stage.inst_exponent]

const Vecs = AbstractVector{<: AbstractVector}

exponents_history(solver::Union{LESolverRecFTLE,
                                AbstractRenormalizer{<: LESolRecFTLE}}) =
    exponents_history(solver.sol)

exponents_history(sol::LESolRecFTLE) = exponents_history(sol.ftle_history)

function exponents_history(ftle_history::Vecs)
    num_attr = length(ftle_history)
    dim_lyap = length(ftle_history[1])

    m = VecMean(dim_lyap)
    s = OnlineStats.Series(m)
    le_hist = similar(ftle_history[1], (dim_lyap, num_attr))
    for i in 1:num_attr
        OnlineStats.fit!(s, ftle_history[i])
        le_hist[:, i] .= mean(m)
    end
    return le_hist
end

exponents_stat_history(sol::LESolRecFTLE) =
    exponents_stat_history(sol.ftle_history)

function exponents_stat_history(ftle_history::Vecs, coverageprob = 0.95)
    num_attr = length(ftle_history)
    dim_lyap = length(ftle_history[1])
    c = cquantile(Normal(), (1 - coverageprob) / 2)
    # Usually num_attr is big so using TDist does not change much here.

    v = VecVariance(dim_lyap)
    s = OnlineStats.Series(v)
    le_hist = similar(ftle_history[1], (dim_lyap, num_attr))
    ci_hist = similar(le_hist)

    OnlineStats.fit!(s, ftle_history[1])
    le_hist[:, 1] .= mean(v)
    ci_hist[:, 1] .= NaN
    for i in 2:num_attr
        OnlineStats.fit!(s, ftle_history[i])
        # c = cquantile(TDist(i - 1), (1 - coverageprob) / 2)
        le_hist[:, i] .= mean(v)
        ci_hist[:, i] .= std(v) / sqrt(i) * c
    end
    return le_hist, ci_hist
end
# MAYBE: use OnlineStats.Bootstrap
# http://joshday.github.io/OnlineStats.jl/latest/api.html#OnlineStats.Bootstrap
