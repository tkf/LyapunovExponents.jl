"""
Auto-generated dynamics for solving phase and tangent dynamics together.
"""
struct PhaseTangentDynamics{F, J}
    "A callable in `(du, u, p, t)` format (inplace version for `ODEProblem`)."
    phase_dynamics::F

    "Temporary variable to store the Jacobian of the `phase_dynamics`."
    jacobian::J

    # TODO: store ForwardDiff config here
end

function fdiff_tangent_dynamics(phase_dynamics, x0)
    dimension = size(x0, 1)
    return PhaseTangentDynamics(
        phase_dynamics,
        similar(x0, (dimension, dimension)),
    )
end

"""Co-evolve phase- and tangent-sapce dynamics."""
@inline function (dyn::PhaseTangentDynamics)(du, u, args...)
    ForwardDiff.jacobian!(
        dyn.jacobian,
        (dx, x) -> dyn.phase_dynamics(dx, x, args...),
        (@view du[:, 1]),       # phase space derivative `dx` goes here
        (@view u[:, 1]),        # current phase space state `x`
    )
    A_mul_B!((@view du[:, 2:end]), dyn.jacobian, (@view u[:, 2:end]))
end
