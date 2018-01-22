"""
The Problem type represents the setup of the Lyapunov exponents calculation.
"""
abstract type AbstractLEProblem end

"""
The Relaxer type represents the calculation required for throwing away
the transient part of the dynamics.
"""
abstract type AbstractRelaxer end

"""
The Solver type represents the core Lyapunov Exponents (LE)
calculation.  The LE calculation is done by calling the in-place
mutating function [`solve!`](@ref).

The methods connecting the three principal types (Problem, Relaxer and
Solver) for the LE calculation are shown in the following diagram:

┌─ Problem ([`AbstractLEProblem`](@ref))                               \\\n
│     │                                                                \\\n
│     │ [`get_relaxer`](@ref), [`relaxed`](@ref)                       \\\n
│     ▼                                                                \\\n
│   Relaxer ([`AbstractRelaxer`](@ref)) ┄┄ ⟲ [`relax!`](@ref)        \\\n
│     │                                                                \\\n
│     │ [`LESolver`](@ref)                                             \\\n
│     │                                                                \\\n
│┄┄┄┄┄┄┄ [`init`](@ref), [`solve`](@ref)                         \\\n
│     │                                                                \\\n
│     ▼                                                                \\\n
└▶ Solver ([`AbstractLESolver`](@ref)) ┄┄ ⟲ [`solve!`](@ref),
                                                 [`step!`](@ref)

"""
abstract type AbstractLESolver end

"""
    LEProblem(phase_prob;  <keyword arguments>)

# Arguments
- `phase_prob`: Phase space dynamics represented in the form of
  `ODEProblem` or `DiscreteProblem` from DifferentialEquations.jl.
  `phase_prob.tspan` represents the inter-orthonormalization-interval.
- `num_tran::Integer`: Number of iterations to through away to get rid
  of the transient dynamics.
- `dim_lyap::Integer`: Number of Lyapunov exponents to be calculated.
  Default to the full system dimension.
- `Q0::Array`: The initial guess of the Gram-Schmidt "Lyapunov vectors".
  Default to the identity matrix.
- `tangent_dynamics::Function`: A vector field for solving phase space
   evolution *and* tangent space evolution together.  If this is not
   provided, `tangent_dynamics` is derived from `phase_prob.f`.  See
   also [`PhaseTangentDynamics`](@ref).
"""
struct LEProblem{DEP} <: AbstractLEProblem
    phase_prob::DEP
    num_tran
    dim_lyap
    Q0
    tangent_dynamics!

    function LEProblem(
            phase_prob::DEP;
            num_tran=1,
            dim_lyap=dimension(phase_prob),
            Q0=eye(dimension(phase_prob), dim_lyap),
            tangent_dynamics=nothing,
            ) where {DEP}
        new{DEP}(phase_prob, num_tran, dim_lyap, Q0, tangent_dynamics)
    end
end

struct Relaxer{LEP} <: AbstractRelaxer
    prob::LEP
    integrator
end

"""
    LESolver(integrator; <keyword arguments>)

A type representing the main calculation of Lyapunov Exponents (LE).
This struct holds all temporary state required for LE calculation.
"""
mutable struct LESolver{Intr} <: AbstractLESolver
    integrator::Intr
    exponents
    num_orth
    phase_state
    tangent_state
    sign_R

    function LESolver(
            integrator::Intr;
            phase_state=get_phase_state(integrator),
            tangent_state=get_tangent_state(integrator),
            dim_lyap=length(phase_state),
            ) where {Intr}
        num_orth = 0
        exponents = zeros(eltype(phase_state), dim_lyap)
        sign_R = Array{Bool}(dim_lyap)
        new{Intr}(
            integrator,
            exponents,
            num_orth,
            phase_state,
            tangent_state,
            sign_R,
        )
    end
end

Base.show(io::IO, solver::LESolver) =
    print(io,
          "#Orth.: ", solver.num_orth, ", ",
          "LEs: ", lyapunov_exponents(solver))
