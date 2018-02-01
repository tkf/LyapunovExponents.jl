module HenonMap
export henon_map
using ..ExampleBase: LEDemo, DiscreteExample

"""
Return a [`LEDemo`](@ref) for the Hénon map.

* M. Hénon, Commun. Math. Phys. Phys. 50, 69-77 (1976)
* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* <https://en.wikipedia.org/wiki/H%C3%A9non_map>
"""
function henon_map(;
        u0=[0.1, 0.1],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=1e-2,
        kwargs...)
    @inline function phase_dynamics!(u_next, u, p, t)
        u_next[1] = 1 + u[2] - 1.4 * u[1]^2
        u_next[2] = 0.3 * u[1]
    end
    @inline @views function tangent_dynamics!(u_next, u, p, t)
        phase_dynamics!(u_next[:, 1], u[:, 1], p, t)
        u_next[1, 2:end] .= u[2, 2:end] .- 1.4 * 2 * u[1, 1] .* u[1, 2:end]
        u_next[2, 2:end] .= 0.3 .* u[1, 2:end]
    end
    LEDemo(DiscreteExample(
        "Hénon map",
        phase_dynamics!, u0, tspan, nothing,
        tangent_dynamics!,
        num_attr,
        [0.41922, -1.62319],   # known_exponents
        atol, rtol,
    ); kwargs...)
end

end
