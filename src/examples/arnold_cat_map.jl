module ArnoldCatMap
export arnold_cat_map
using ..ExampleBase: LEDemo, DiscreteExample

@inline function phase_dynamics!(u_next, u, M, t)
    A_mul_B!(u_next, M, u)
    u_next .= mod.(u_next, 1)
end

@inline @views function tangent_dynamics!(u_next, u, M, t)
    phase_dynamics!(u_next[:, 1], u[:, 1], M, t)
    for i in 2:size(u)[2]
        @inbounds A_mul_B!(u_next[:, i], M, u[:, i])
    end
end

"""
Return a [`LEDemo`](@ref) for the Arnold's cat map

* <https://en.wikipedia.org/wiki/Arnold%27s_cat_map>
"""
function arnold_cat_map(;
        u0::AbstractVector = [0.1, 0.1],
        M::AbstractMatrix = [2 1; 1 1],
        tspan=10,
        t_attr=10000,
        atol=0, rtol=1e-4,
        kwargs...)
    @assert (length(u0), length(u0)) == size(M)
    LEDemo(DiscreteExample(
        "Arnold's cat map",
        phase_dynamics!, u0, tspan, M,
        tangent_dynamics!,
        t_attr,
        log.(sort!(eigvals(M), rev=true)), # known_exponents
        atol, rtol,
    ); kwargs...)
end

end
