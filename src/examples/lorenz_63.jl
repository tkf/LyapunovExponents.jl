module Lorenz63
export lorenz_63
using ..ExampleBase: LEDemo, ContinuousExample

"""
Return a [`LEDemo`](@ref) for the Lorenz system.

* <https://en.wikipedia.org/wiki/Lorenz_system>
* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* E. N. Lorenz, J. Atmos. Sci. 20, 130-141 (1963)
"""
function lorenz_63(;
        u0=[0.1, 0.1, 0.1],
        tspan=(0.0, 1.0),
        num_attr=4000,
        atol=0, rtol=1e-2,
        kwargs...)
    @inline function phase_dynamics!(du, u, p, t)
        du[1] = 10.0(u[2]-u[1])
        du[2] = u[1]*(28.0-u[3]) - u[2]
        du[3] = u[1]*u[2] - (8/3)*u[3]
    end
    @inline @views function tangent_dynamics!(du, u, p, t)
        phase_dynamics!(du[:, 1], u[:, 1], p, t)
        du[1, 2:end] .= 10.0 .* (u[2, 2:end] .- u[1, 2:end])
        du[2, 2:end] .=
            u[1, 2:end] .* (28.0 .- u[3, 1]) .-
            u[1, 1] .* u[3, 2:end] .-
            u[2, 2:end]
        du[3, 2:end] .=
            u[1, 2:end] .* u[2, 1] .+ u[1, 1] .* u[2, 2:end] .-
            (8/3) .* u[3, 2:end]
    end
    LEDemo(ContinuousExample(
        "Lorenz (1963)",
        phase_dynamics!,
        u0, tspan, num_attr,
        tangent_dynamics!,
        [0.9056, 0, -14.5723],  # known_exponents
        atol, rtol,
    ); kwargs...)
end

end
