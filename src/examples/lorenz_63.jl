module Lorenz63
export lorenz_63

using Parameters: @with_kw, @unpack
using ..ExampleBase: LEDemo, ContinuousExample

@with_kw struct Lorenz63Param
    σ::Float64 = 10.0
    ρ::Float64 = 28.0
    β::Float64 = 8/3
end

function phase_dynamics!(du, u, p, t)
    @unpack σ, ρ, β = p
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β*u[3]
end

@inline @views function tangent_dynamics!(du, u, p, t)
    @unpack σ, ρ, β = p
    phase_dynamics!(du[:, 1], u[:, 1], p, t)
    du[1, 2:end] .= σ .* (u[2, 2:end] .- u[1, 2:end])
    du[2, 2:end] .=
        u[1, 2:end] .* (ρ .- u[3, 1]) .-
        u[1, 1] .* u[3, 2:end] .-
        u[2, 2:end]
    du[3, 2:end] .=
        u[1, 2:end] .* u[2, 1] .+ u[1, 1] .* u[2, 2:end] .-
        β .* u[3, 2:end]
end

"""
Return a [`LEDemo`](@ref) for the Lorenz system.

* <https://en.wikipedia.org/wiki/Lorenz_system>
* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* E. N. Lorenz, J. Atmos. Sci. 20, 130-141 (1963)
"""
function lorenz_63(;
        u0=[0.1, 0.1, 0.1],
        t_renorm=1.0,
        t_attr=100000,
        atol=0, rtol=1e-2,
        # terminator_options = [:atol => atol, :rtol => rtol],  # TODO: use it
        terminator_options = [],
        kwargs...)
    LEDemo(ContinuousExample(
        "Lorenz (1963)",
        phase_dynamics!, u0, t_renorm, Lorenz63Param(),
        tangent_dynamics!,
        t_attr,
        [0.9056, 0, -14.5723],  # known_exponents
        atol, rtol,
        terminator_options,
    ); kwargs...)
end

end
