"""
Auto-generated dynamics for solving phase and tangent dynamics together.
"""
type PhaseTangentDynamics
    "A callable in `(t, u, du)` format (inplace version for `ODEProblem`)."
    phase_dynamics!

    "Temporary variable to store the Jacobian of the `phase_dynamics!`."
    jacobian

    # TODO: store ForwardDiff config here

    function PhaseTangentDynamics(phase_dynamics!, array_type,
                                  dimension::Integer)
        jacobian = array_type(dimension, dimension)
        new(phase_dynamics!,
            jacobian,
            )
    end
end

function PhaseTangentDynamics(
        phase_dynamics!, x0::AbstractArray{T, 1}) where {T}
    PhaseTangentDynamics(phase_dynamics!, Array{eltype(x0)}, length(x0))
end

function PhaseTangentDynamics(
        phase_dynamics!, u0::AbstractArray{T, 2}) where {T}
    PhaseTangentDynamics(phase_dynamics!, u0[:, 1])
end

"""Co-evolve phase- and tangent-sapce dynamics."""
@inline function (dyn::PhaseTangentDynamics)(t, u, du)
    ForwardDiff.jacobian!(
        dyn.jacobian,
        (dx, x) -> dyn.phase_dynamics!(t, x, dx),
        (@view du[:, 1]),       # phase space derivative `dx` goes here
        (@view u[:, 1]),        # current phase space state `x`
    )
    A_mul_B!((@view du[:, 2:end]), dyn.jacobian, (@view u[:, 2:end]))
end
