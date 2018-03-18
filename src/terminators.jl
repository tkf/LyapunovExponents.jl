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
    cov::Float64
end
struct StableConvDetail <: ConvDetail
    high::Float64
    low::Float64
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
                       i::Int, err, th, var, cov)
    push!(history.errors[i], err)
    push!(history.thresholds[i], th)
    push!(history.details[i], UnstableConvDetail(var, cov))
end

function record_error!(history::ConvergenceHistory,
                       ::Type{Val{StableConvError}},
                       i::Int, err, th, high, low)
    push!(history.errors[i], err)
    push!(history.thresholds[i], th)
    push!(history.details[i], StableConvDetail(high, low))
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
    max_cutoff_corr::Float64
    next_check::Int
end
# TODO: separate state (next_check) and settings.

function AutoCovTerminator(;
        rtol = 0.1,
        atol = 0.1,
        chunks = 2,
        max_mle = 0.01,
        max_cutoff_corr = 0.05,
        first_check = 100,
        )
    @assert chunks > 1
    AutoCovTerminator(rtol, atol,
                      chunks, max_mle,
                      max_cutoff_corr,
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
        ok = ok && err < th
        ok = ok && abs(c[end]) < abs(c[1]) * tmnr.max_cutoff_corr
        record_error!(history, Val{UnstableConvError}, i, err, th, c[1], c[end])
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
    ftle = @view sol.ftle_history[1:n]
    dim_lyap = size(ftle[1], 1)

    history = sol.convergence
    new_error!(history, n, StableConvError)
    rerr = 1.0
    ok = true
    for i in 1:dim_lyap
        w = floor(Int, n / tmnr.chunks)
        ms = [mean(x[i] for x in ftle[j:j+w-1]) for j in 1:w:n-w+1]
        @assert length(ms) == tmnr.chunks
        h = maximum(ms)
        l = minimum(ms)
        th = atol + max(abs(h), abs(l)) * rtol
        err = h - l
        rerr = max(rerr, err / th)
        ok = ok && err < th
        record_error!(history, Val{StableConvError}, i, err, th, h, l)
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


"""Choose appropriate terminator."""
function Terminator(::Void, prob, sol; tmnr_opts...)
    if sol isa LESolution{true, true}
        tmnr = AutoCovTerminator(; tmnr_opts...)
    else
        @assert isempty(tmnr_opts)
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
