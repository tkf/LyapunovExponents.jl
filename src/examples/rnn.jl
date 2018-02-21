module RNN
export beer_95

using Parameters: @unpack
using ..ExampleBase: LEDemo, ContinuousExample

σ(x) = 1 / (1 + exp(-x))
σ′(x) = σ(x) * (1 - σ(x))

type ContinuousRNN
    w
    θ
    τ
end

@inline function phase_dynamics!(du, u, p, t)
    @unpack w, θ, τ = p
    du .= (- u .+ w * σ.(u .+ θ)) ./ τ
end

@inline @views function tangent_dynamics!(du, u, p, t)
    @unpack w, θ, τ = p
    phase_dynamics!(du[:, 1], u[:, 1], p, t)
    Y = u[:, 2:end]
    du[:, 2:end] .= (- Y .+ w * Diagonal(σ′.(u[:, 1] .+ θ)) * Y) ./ τ
end

"""
Return a [`LEDemo`](@ref) for a low-dimensional chaotic
continuous-time recurrent neural networks by Beer (1995).

* Beer, R. D. (1995). On the dynamics of small continuous-time recurrent
  neural networks. Adapt. Behav., 3(4), 469–509.
  <https://doi.org/10.1177/105971239500300405>.
  (Figure 9D)
"""
function beer_95(;
        u0=[0.1, 0.1, 0.1],
        tspan=10.0,
        t_attr=4000.0,
        atol=0, rtol=5e-1,
        kwargs...)
    param = ContinuousRNN(
        # w
        [  5.422  -0.018  2.75
          -0.24    4.59   1.21
           0.535  -2.25   3.885 ],
        # θ
        [-4.108, -2.787, -1.114],
        # τ
        [1.0, 2.5, 1.0],
    )
    LEDemo(ContinuousExample(
        "Beer (1995)",
        phase_dynamics!, u0, tspan, param,
        tangent_dynamics!,
        t_attr,
        [0.010, 0],   # known_exponents
        atol, rtol,
    ); kwargs...)
end

end
