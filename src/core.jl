using DiffEqBase: DEProblem, remake
using Distributions: Normal, cquantile

dimension(prob::DEProblem) = length(prob.u0)

PhaseRelaxer(prob::LEProblem,
             ::LEProblem,
             ::LESolution,
             ::Terminator,
             integrator_options) =
    PhaseRelaxer(get_integrator(prob.phase_prob; integrator_options...),
                 prob.t_tran)

function step!(relaxer::PhaseRelaxer)
    step!(relaxer.integrator, relaxer.t_tran)
    relaxer.i += 1
end
Base.length(relaxer::PhaseRelaxer) = 1  # FIXME

function get_tangent_prob(prob;
                          x0 = prob.phase_prob.u0,
                          kwargs...)
    u0 = copy(prob.tangent_prob.u0)
    u0[:, 1] = x0
    return remake(prob.tangent_prob; u0=u0, kwargs...)
end

current_state(relaxer::Union{PhaseRelaxer, AbstractRenormalizer}) =
    current_state(relaxer.integrator)

function get_tangent_integrator(prob::LEProblem, relaxer;
                                integrator_options...)
    # Carry on simulated time to tangent integrator [*]:
    tspan = if relaxer.integrator isa ODEIntegrator
        (relaxer.integrator.t, Inf)
    else
        # Support this in discrete problem
        prob.phase_prob.tspan
    end
    x0 = current_state(relaxer)
    tangent_prob = get_tangent_prob(prob; x0=x0, tspan=tspan)
    return get_integrator(tangent_prob; integrator_options...)
end
# [*] Another way to propagate simulated time to tangent integrator is
# to use set_t!() after init().  However, init() calls phase_prob.f
# with phase_prob.tspan[1] which is less than relaxer.integrator.t.
# It causes problems in some rare cases with stateful problem
# definitions.

TangentRenormalizer(relaxer::PhaseRelaxer, prob::LEProblem, sol::LESolution,
                    tmnr::Terminator, integrator_options) =
    TangentRenormalizer(get_tangent_integrator(prob, relaxer;
                                               integrator_options...),
                        prob.t_attr, prob.t_renorm, sol, tmnr)

MLERenormalizer(relaxer::PhaseRelaxer, prob::LEProblem, sol::LESolution,
                tmnr::Terminator, integrator_options) =
    MLERenormalizer(get_tangent_integrator(prob, relaxer;
                                           integrator_options...),
                    prob.t_attr, prob.t_renorm, sol, tmnr)

Base.length(stage::AbstractRenormalizer) =
    ceil(Int, stage.t_attr / stage.t_renorm)

@inline function keepgoing!(stage::AbstractRenormalizer)
    u0 = current_state(stage)
    u0[:, 1] = stage.phase_state
    u0[:, 2:end] = stage.tangent_state
    set_u!(stage.integrator, u0)
    step!(stage.integrator, stage.t_renorm, true)
    assert_success(stage.integrator)
    # TODO: maybe use step!(solver.integrator, dt) aka "inexact step"
    u0 = current_state(stage)
    stage.phase_state[:] = @view u0[:, 1]
    stage.tangent_state[:] = @view u0[:, 2:end]
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
    step!(stage::AbstractRenormalizer)

Evolve the dynamics and then do an orthonormalization.
"""
function step!(stage::AbstractRenormalizer)
    keepgoing!(stage)
    post_evolve!(stage)
    record!(stage, Val{:OS})
    record!(stage, Val{:FTLE})
    stage.is_finished = should_terminate!(stage.tmnr, stage.sol)
end

is_finished(stage::AbstractRenormalizer) =
    stage.is_finished || (stage.i >= length(stage))

function post_evolve!(stage::TangentRenormalizer)
    dim_lyap = length(stage.inst_exponents)
    P = stage.tangent_state
    F = qrfact!(P)
    Q = stage.Q
    A_mul_B!(F[:Q], eye!(Q))  # Q = Matrix(F[:Q])[...]; but lesser allocation
    R = stage.R = F[:R]

    sign_R = stage.sign_R
    sign_diag!(sign_R, R)        # sign_R = diagm(sign(diag(R)))
    A_mul_sign!(Q, sign_R)       # Q = Q * sign_R
    sign_mul_A!(sign_R, R)       # R = signR * R

    stage.i += 1
    stage.tangent_state, stage.Q = Q, stage.tangent_state
    # At this point:
    # stage.tangent_state === Q == 𝑮ₙ₊ₖ
    # stage.Q is a mess, as qrfact! use it as a "buffer"
end

function post_evolve!(stage::MLERenormalizer)
    v = stage.tangent_state[:, 1]
    r = norm(v)
    n = (stage.i += 1)
    stage.exponent = lyap_add_r(n, stage.exponent, r)
    stage.inst_exponent = log(r) / t_chunk(stage)
    # TODO: don't re-calculate log(r)
    stage.tangent_state[:, 1] .= v ./ r
end

t_chunk(stage) = stage.t_renorm

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

function record!(stage::AbstractRenormalizer{<: LESolRecFTLE},
                 ::Type{Val{:FTLE}})
    stage.sol.ftle_history[stage.i] .= ftle(stage)
    stage.sol.num_orth = stage.i
end

current_result(stage::TangentRenormalizer) = stage.inst_exponents
current_result(stage::MLERenormalizer) = stage.inst_exponent

"""
    lyapunov_exponents(solver)

Get the result of Lyapunov exponents calculation stored in `solver`.
"""
lyapunov_exponents(solver::Union{LESolverRecOS,
                                 AbstractRenormalizer{<: LESolRecOS}}) =
    lyapunov_exponents(solver.sol)

lyapunov_exponents(sol::LESolRecOS) = mean(sol.main_stat)

@inline function lyapunov_exponents(stage::MLERenormalizer)
    [stage.exponent / t_chunk(stage)]
end

"""Get finite-time Lyapunov exponents (FTLE)"""
ftle(stage::TangentRenormalizer) = stage.inst_exponents
ftle(stage::MLERenormalizer) = [stage.inst_exponent]

ftle_history(solver::Union{LESolverRecFTLE,
                           AbstractRenormalizer{<: LESolRecFTLE}},
             args...) =
    ftle_history(solver.sol, args...)

ftle_history(sol::LESolRecFTLE) = @view sol.ftle_history[1:sol.num_orth]
ftle_history(sol::LESolRecFTLE, i) =
    [x[i] for x in sol.ftle_history[1:sol.num_orth]]

function ftle_history(::Type{A}, sol::LESolRecFTLE) where {A <: AbstractMatrix}
    dim_lyap = length(sol.ftle_history[1])
    num_orth = sol.num_orth
    arr = A(dim_lyap, num_orth)
    for (i, ftle) in enumerate(ftle_history(sol))
        arr[:, i] .= ftle
    end
    return arr
end

const Vecs = AbstractVector{<: AbstractVector}

exponents_history(sol) = exponents_history(ftle_history(sol))

function exponents_history(ftle_history::Vecs)
    t_attr = length(ftle_history)
    dim_lyap = length(ftle_history[1])

    m = VecMean(dim_lyap)
    s = OnlineStats.Series(m)
    le_hist = similar(ftle_history[1], (dim_lyap, t_attr))
    for i in 1:t_attr
        OnlineStats.fit!(s, ftle_history[i])
        le_hist[:, i] .= mean(m)
    end
    return le_hist
end

exponents_stat_history(sol::Union{LESolRecFTLE,
                                  AbstractRenormalizer{<: LESolRecFTLE},
                                  LESolverRecFTLE}) =
    exponents_stat_history(ftle_history(sol))

function exponents_stat_history(ftle_history::Vecs, coverageprob = 0.95)
    t_attr = length(ftle_history)
    dim_lyap = length(ftle_history[1])
    c = cquantile(Normal(), (1 - coverageprob) / 2)
    # Usually t_attr is big so using TDist does not change much here.

    v = VecVariance(dim_lyap)
    s = OnlineStats.Series(v)
    le_hist = similar(ftle_history[1], (dim_lyap, t_attr))
    ci_hist = similar(le_hist)

    OnlineStats.fit!(s, ftle_history[1])
    le_hist[:, 1] .= mean(v)
    ci_hist[:, 1] .= NaN
    for i in 2:t_attr
        OnlineStats.fit!(s, ftle_history[i])
        # c = cquantile(TDist(i - 1), (1 - coverageprob) / 2)
        le_hist[:, i] .= mean(v)
        ci_hist[:, i] .= std(v) / sqrt(i) * c
    end
    return le_hist, ci_hist
end
# MAYBE: use OnlineStats.Bootstrap
# http://joshday.github.io/OnlineStats.jl/latest/api.html#OnlineStats.Bootstrap
