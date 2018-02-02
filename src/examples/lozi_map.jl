module LoziMap
export lozi_map

using Parameters: @with_kw, @unpack
using ..ExampleBase: LEDemo, DiscreteExample

@with_kw struct LoziMapParam
    a::Float64 = 1.4
    b::Float64 = 0.3
end

@inline function phase_dynamics!(u_next, u, p, t)
    @unpack a, b = p
    u_next[1] = 1 + u[2] - a * abs(u[1])
    u_next[2] = b * u[1]
end

@inline @views function tangent_dynamics!(u_next, u, p, t)
    @unpack a, b = p
    phase_dynamics!(u_next[:, 1], u[:, 1], p, t)
    u_next[1, 2:end] .= u[2, 2:end] .- a * sign(u[1, 1]) .* u[1, 2:end]
    u_next[2, 2:end] .= b .* u[1, 2:end]
end

"""
Return a [`LEDemo`](@ref) for the Lozi map.
"""
function lozi_map(;
        u0=[0.1, 0.1],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=0.2,
        kwargs...)
    LEDemo(DiscreteExample(
        "Lozi map",
        phase_dynamics!, u0, tspan, LoziMapParam(),
        tangent_dynamics!,
        num_attr,
        Float64[],   # known_exponents
        atol, rtol,
    ); kwargs...)
end

end