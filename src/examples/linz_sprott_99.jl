module LinzSprott99
export linz_sprott_99
using OrdinaryDiffEq
using ..ExampleBase: LEDemo, ContinuousExample

@inline function phase_dynamics!(du, u, A, t)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -A * u[3] - u[2] - (u[1] > 0 ? 1 : -1) * u[1] + 1
end

@inline @views function tangent_dynamics!(du, u, A, t)
    phase_dynamics!(du[:, 1], u[:, 1], A, t)
    du[1, 2:end] = u[2, 2:end]
    du[2, 2:end] = u[3, 2:end]
    du[3, 2:end] .=
        -A .* u[3, 2:end] .- u[2, 2:end] .-
        (u[1, 1] > 0 ? 1 : -1) .* u[1, 2:end]
end

"""
Return a [`LEDemo`](@ref) for the simplest piecewise linear
dissipative chaotic flow.

* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* S. J. Linz and J. C. Sprott, Phys. Lett. A 259, 240-245 (1999)
"""
function linz_sprott_99(;
        u0=[0.1, 0.1, 0.1],
        t_renorm=10.0,
        t_attr=1000000,
        atol=1e-2, rtol=1e-2,
        terminator_options = [:atol => atol, :rtol => rtol],
        integrator_options = [
            :alg => Vern6(),
            # :alg => BS5(),
        ],
        kwargs...)
    A = 0.6
    LEDemo(ContinuousExample(
        "Linz & Sprott (1999)",
        phase_dynamics!, u0, t_renorm, A,
        tangent_dynamics!,
        t_attr,
        [0.0362, 0, -0.6362],   # known_exponents
        atol, rtol,
        terminator_options,
        integrator_options = integrator_options,
    ); kwargs...)
end

end
