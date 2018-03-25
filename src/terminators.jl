using StatsBase: autocov

"""
    check_cache_type(terminator, prob, sol)
"""
check_cache_type(::T, ::P, ::S) where {T, P, S} =
    error("Terminator type is incompatible with the ",
          "problem and solution types.\n",
          "terminator: $T\n",
          "problem: $P\n",
          "solution: $S")

should_terminate!(::T, ::S) where {T, S} =
    error("`should_terminate!` not defined for $T and $S")


struct NullTerminator <: Terminator end
check_cache_type(::NullTerminator, ::Any, ::LESolution) = true
should_terminate!(::NullTerminator, ::LESolution) = false

@enum ConvErrorType UnstableConvError=1 StableConvError=2

abstract type ConvDetail end
struct UnstableConvDetail <: ConvDetail
    var::Float64
    tail_cov::Float64
    tail_ok::Bool
end

struct ConvergenceHistory
    orth::Vector{Int}
    kinds::Vector{ConvErrorType}
    errors::Vector{Vector{Float64}}
    thresholds::Vector{Vector{Float64}}
    details::Vector{Vector{ConvDetail}}
end

ConvergenceHistory(dim_lyap::Int) =
    ConvergenceHistory(
        Vector{Int}(),
        Vector{ConvErrorType}(),
        [Vector{Float64}() for _ in 1:dim_lyap],
        [Vector{Float64}() for _ in 1:dim_lyap],
        [Vector{ConvDetail}() for _ in 1:dim_lyap],
    )

function new_error!(history::ConvergenceHistory, n::Int, kind::ConvErrorType)
    push!(history.orth, n)
    push!(history.kinds, kind)
end

function record_error!(history::ConvergenceHistory,
                       ::Type{Val{UnstableConvError}},
                       i::Int, err, th, args...)
    push!(history.errors[i], err)
    push!(history.thresholds[i], th)
    push!(history.details[i], UnstableConvDetail(args...))
end

function record_error!(history::ConvergenceHistory,
                       i::Int, err, th, detail)
    push!(history.errors[i], err)
    push!(history.thresholds[i], th)
    push!(history.details[i], detail)
end

function assert_consistent(history::ConvergenceHistory)
    @assert [length(history.orth)] ==
        [length(history.kinds)] ==
        unique(length.(history.errors)) ==
        unique(length.(history.thresholds)) ==
        unique(length.(history.details))
end

absolute_error(history::ConvergenceHistory, i::Int) = history.errors[i]

relative_error(history::ConvergenceHistory, i::Int) =
    history.errors[i] ./ history.thresholds[i]


mutable struct AutoCovTerminator <: Terminator
    rtol::Float64
    atol::Float64
    chunks::Int
    max_mle::Float64
    max_tail_corr::Float64
    tail_ratio::Float64
    max_fp_cv::Float64
    min_periodic_corr::Float64
    skip_ratio::Float64
    next_check::Int
end
# TODO: separate state (next_check) and settings.
#
# TODO: Separate terminator for unstable (chaotic) and stable cases
# and add a terminator combinator for composing them.  The stable case
# can be decomposed into more sub-components (see stable_le_error).
#
# TODO: rtol and atol are used for unstable (chaotic) and stable cases
# but they have different semantics in each case.  Find a way to
# translate one to the other and/or use different option.  Using
# different option would be trivial once the combinator is
# implemented.

function AutoCovTerminator(;
        rtol = 0.1,
        atol = 0.1,
        chunks = 2,
        max_mle = 0.01,
        max_tail_corr = 0.05,
        tail_ratio = 0.1,
        skip_ratio = 0.1,
        max_fp_cv = 0.001,
        min_periodic_corr = 0.8,
        first_check = 400,
        )
    @assert chunks > 1
    AutoCovTerminator(rtol, atol,
                      chunks, max_mle,
                      max_tail_corr,
                      tail_ratio,
                      skip_ratio,
                      max_fp_cv,
                      min_periodic_corr,
                      first_check)
end

check_cache_type(::AutoCovTerminator, ::Any, ::LESolution{true, true}) = true
# Need RecOS=true and RecFTLE=true.


function should_terminate!(tmnr::AutoCovTerminator, sol::LESolRecFTLE)
    atol = tmnr.atol
    rtol = tmnr.rtol
    n = sol.num_orth
    if n < tmnr.next_check
        return false
    end

    les = lyapunov_exponents(sol)
    if les[1] < tmnr.max_mle
        return should_terminate_stable!(tmnr, sol)
    end

    ftle = @view sol.ftle_history[1:n]
    dim_lyap = size(ftle[1], 1)

    history = sol.convergence
    new_error!(history, n, UnstableConvError)
    rerr = 1.0
    ok = true
    for i in 1:dim_lyap
        xs = [x[i] for x in ftle]
        m = mean(xs)
        err, c = correlated_sem_with_cov(xs)
        th = atol + abs(m) * rtol
        rerr = max(rerr, err / th)
        tail_width = ceil(Int, length(c) * tmnr.tail_ratio)
        tail_cov = mean(abs, @view c[end-tail_width:end])
        tail_ok = tail_cov < abs(c[1]) * tmnr.max_tail_corr
        ok = ok && err < th
        ok = ok && tail_ok
        record_error!(history, Val{UnstableConvError}, i, err, th,
                      c[1], tail_cov, tail_ok)
    end
    assert_consistent(history)

    if !ok
        # aiming for err < th/2 (so that's why "4 *")
        tmnr.next_check = ceil(Int, 4 * rerr^2 * n)
        tmnr.next_check = min(tmnr.next_check, length(sol.ftle_history))
    end
    sol.converged = ok
    return ok
end


function should_terminate_stable!(tmnr::AutoCovTerminator, sol::LESolRecFTLE)
    atol = tmnr.atol
    rtol = tmnr.rtol
    n = sol.num_orth
    skip = max(1, floor(Int, n * tmnr.skip_ratio))
    dim_lyap = size(sol.ftle_history[1], 1)

    history = sol.convergence
    new_error!(history, n, StableConvError)
    rerr = 1.0
    ok = true
    for i in 1:dim_lyap
        ftle = @view ftle_history(sol, i)[skip:end]
        m = mean(ftle)
        th = atol + abs(m) * rtol
        err, detail = stable_le_error(ftle)
        if isfinite(err)
            rerr = max(rerr, err / th)
        end
        ok = ok && err < th
        record_error!(history, i, err, th, detail)
    end
    assert_consistent(history)

    if !ok
        # aiming for err < th/2 (so that's why "2 *")
        tmnr.next_check = ceil(Int, 2 * rerr * n)
        tmnr.next_check = min(tmnr.next_check, length(sol.ftle_history))
    end
    sol.converged = ok
    return ok
end


"""
    correlated_sem_with_cov(xs [, cutoff ≈ √length(xs)])

Standard error of mean of (possibly) correlated data.
https://stats.stackexchange.com/a/274672
"""
function correlated_sem_with_cov(xs::AbstractVector,
                                 cutoff = ceil(Int, sqrt(length(xs))))
    c = autocov(xs, 0:cutoff)
    σ = c[1]
    g = @view c[2:end]
    n = length(xs)
    k = 1:cutoff
    sem² = (σ + 2sum((n - k) ./ n .* g)) / n
    sem = sqrt(max(σ / n, sem²))    # TODO: should I use max?
    return (sem, c)
end


struct FixedPointConvDetail <: ConvDetail
    var::Float64
end

struct NonNegativeAutoCovConvDetail <: ConvDetail
    var::Float64
    min_corr::Float64
end

struct NonPeriodicConvDetail <: ConvDetail
    var::Float64
    max_corr::Float64
end

struct PeriodicConvDetail <: ConvDetail
    period::Int
end

stable_le_error(ftle, opts) =
    stable_le_error(ftle;
                    max_fp_cv = opts.max_fp_cv,
                    min_periodic_corr = opts.min_periodic_corr,
                    )

function stable_le_error(ftle::AbstractVector;
                         max_fp_cv = 0.001,
                         min_periodic_corr = 0.8,
                         cutoff = length(ftle) ÷ 2)

    v = var(ftle, corrected=false)  # = autocov(ftle, 0:0)[1]
    m = mean(ftle)

    # Handle non-fluctuating LE first.  Likely a fixed point.
    if sqrt(v) / abs(m) < max_fp_cv
        # TODO: detect fixed-point directly in the phase space.
        return 0.0, FixedPointConvDetail(v)
    end

    # Detect "period" in ftle.  Here, "period" is refered to the index
    # with maximum covariance *after* it crosses the negative value.
    g = autocov(ftle, 1:cutoff)
    i = findfirst(x -> x < 0, g)
    if i == 0
        # Dynamics assumed to be periodic but autocovariance does not
        # cross 0.  Probably length(ftle) is too small.
        return (Inf, NonNegativeAutoCovConvDetail(v, minimum(g) / v))
    end
    min_cov = min_periodic_corr * v
    j = findfirst(x -> x > min_cov, @view g[i + 1:end])
    if j == 0
        max_cov = maximum(@view g[i + 1:end])
        return (Inf, NonPeriodicConvDetail(v, max_cov / v))
    end
    period = i + j

    n = length(ftle)
    peak = maximum(abs(x - m) for x in ftle)
    err = period / n * peak
    return err, PeriodicConvDetail(period)
end


"""Choose appropriate terminator."""
function Terminator(::Void, prob, sol; terminator_options...)
    return Terminator(sol isa LESolution{true, true},
                      prob, sol; terminator_options...)
end

function Terminator(flag::Bool, prob, sol; terminator_options...)
    if flag
        tmnr = AutoCovTerminator(; terminator_options...)
    else
        @assert isempty(terminator_options)
        tmnr = NullTerminator()
    end
    check_cache_type(tmnr, prob, sol)
    return tmnr
end

function Terminator(tmnr::Terminator, prob, sol)
    check_cache_type(tmnr, prob, sol)
    return tmnr
end


has_convergence_history(::Any, i::Int) = false
has_convergence_history(solver::LESolver, i) =
    has_convergence_history(solver.sol, i)
has_convergence_history(sol::LESolution, i) =
    has_convergence_history(sol.convergence, i)
has_convergence_history(history::ConvergenceHistory, i::Int) =
    1 <= i <= length(history.errors)

convergence_history(solver::LESolver) = convergence_history(solver.sol)
convergence_history(sol::LESolution) = sol.convergence
