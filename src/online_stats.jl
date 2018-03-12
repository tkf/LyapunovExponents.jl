using OnlineStats: ExactStat, VectorOb, smooth!, smooth, value, unbias
import OnlineStats
import OnlineStatsBase

struct VecMean <: ExactStat{1}
    μs::Vector{Float64}
    VecMean(μs::Vector{Float64}) = new(μs)
end
VecMean(len::Integer) = VecMean(zeros(Float64, len))
OnlineStatsBase.fit!(o::VecMean, x::VectorOb, γ::Float64) = smooth!(o.μs, x, γ)
Base.mean(o::VecMean) = value(o)


mutable struct VecVariance <: ExactStat{1}
    σ2s::Vector{Float64}     # biased variance
    μs::Vector{Float64}
    tmp::Vector{Float64}
    nobs::Int
    function VecVariance(σ2s::Vector{Float64}, μs::Vector{Float64})
        @assert size(σ2s) == size(μs)
        new(σ2s, μs, similar(μs), 0)
    end
end
VecVariance(len::Integer) = VecVariance(zeros(Float64, len),
                                        zeros(Float64, len))

function OnlineStatsBase.fit!(o::VecVariance, x::VectorOb, γ::Float64)
    o.nobs += 1
    o.tmp .= smooth.(o.μs, x, γ)
    o.σ2s .= smooth.(o.σ2s, (x .- o.μs) .* (x .- o.tmp), γ)
    o.tmp, o.μs = o.μs, o.tmp
end

OnlineStatsBase.value(o::VecVariance) =
    o.nobs < 2 ? fill(NaN, size(o.μs)) : o.σ2s .* unbias(o)
Base.var(o::VecVariance) = value(o)
Base.std(o::VecVariance) = sqrt.(var(o))
Base.mean(o::VecVariance) = o.μs

# This is required for OnlineStats 0.15.3:
OnlineStats.nobs(o::VecVariance) = o.nobs
