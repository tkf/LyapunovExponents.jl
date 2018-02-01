module ArnoldCatMap
export arnold_cat_map
using ..ExampleBase: LEDemo, DiscreteExample

"""
Return a [`LEDemo`](@ref) for the Arnold's cat map

* <https://en.wikipedia.org/wiki/Arnold%27s_cat_map>
"""
function arnold_cat_map(;
        u0=[0.1, 0.1],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=1e-4,
        kwargs...)
    @assert size(u0) == (2,)
    M = [2 1; 1 1]
    @inline function phase_dynamics!(u_next, u, p, t)
        A_mul_B!(u_next, M, u)
        u_next .= mod.(u_next, 1)
    end
    @inline @views function tangent_dynamics!(u_next, u, p, t)
        phase_dynamics!(u_next[:, 1], u[:, 1], p, t)
        for i in 2:size(u)[2]
            @inbounds A_mul_B!(u_next[:, i], M, u[:, i])
        end
    end
    LEDemo(DiscreteExample(
        "Arnold's cat map",
        phase_dynamics!, u0, tspan, nothing,
        tangent_dynamics!,
        num_attr,
        log.(sort!(eigvals(M), rev=true)), # known_exponents
        atol, rtol,
    ); kwargs...)
end

end
