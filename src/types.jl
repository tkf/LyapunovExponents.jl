using DiffEqBase: ODEProblem, DiscreteProblem

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

<code class=LE-diagram>  \\\n
┌─ Problem ([`AbstractLEProblem`](@ref))                               \\\n
│     │                                                                \\\n
│     │ [`get_relaxer`](@ref), [`relaxed`](@ref)                       \\\n
│     ▼                                                                \\\n
│   Relaxer ([`AbstractRelaxer`](@ref)) ┄┄ ⟲ [`relax!`](@ref)        \\\n
│     │                                                                \\\n
│     │ [`init`](@ref)                                                 \\\n
│     │                                                                \\\n
│┄┄┄┄┄┄┄ [`init`](@ref), [`solve`](@ref)                         \\\n
│     │                                                                \\\n
│     ▼                                                                \\\n
└▶ Solver ([`AbstractLESolver`](@ref)) ┄┄ ⟲ [`solve!`](@ref),
                                                 [`step!`](@ref)         \\\n
</code>

"""
abstract type AbstractLESolver{Intr} end

"""
    LEProblem(phase_prob, num_attr; <keyword arguments>)
    LEProblem(phase_prob; num_attr, <keyword arguments>)

# Arguments
- `phase_prob`: Phase space dynamics represented in the form of
  `ODEProblem` or `DiscreteProblem` from DifferentialEquations.jl.
  `phase_prob.tspan` represents the inter-orthonormalization-interval.
- `num_attr::Integer`: Number of orthonormalizations for calculating
  Lyapunov Exponents.  The simulated time of the system for this
  calculation is given by `num_attr * (tspan[1] - tspan[0])`.
  This argument is always required and can be given as positional or
  keyword argument.
- `num_tran::Integer`: Number of iterations to through away to get rid
  of the transient dynamics.
- `dim_lyap::Integer`: Number of Lyapunov exponents to be calculated.
  Default to the full system dimension.
- `Q0::Array`: The initial guess of the Gram-Schmidt "Lyapunov vectors".
- `tangent_dynamics::Function`: A vector field for solving phase space
   evolution *and* tangent space evolution together.  If this is not
   provided, `tangent_dynamics` is derived from `phase_prob.f`.  See
   also [`PhaseTangentDynamics`](@ref).
"""
struct LEProblem{DEP} <: AbstractLEProblem
    phase_prob::DEP
    num_tran
    num_attr
    dim_lyap
    Q0
    tangent_dynamics!

    function LEProblem{DEP}(
        phase_prob::DEP, num_attr::Int;
            num_tran=1,
            dim_lyap=dimension(phase_prob),
            Q0 = default_Q0(phase_prob, dimension(phase_prob), dim_lyap),
            tangent_dynamics! = nothing,
            ) where {DEP}
        if ! is_semi_unitary(Q0)
            error("Columns in Q0 are not orthonormal.")
        end
        new{DEP}(phase_prob, num_tran, num_attr,
                 dim_lyap, Q0, tangent_dynamics!)
    end
end

LEProblem(phase_prob::DEP, num_attr; kwargs...) where {DEP <: ODEProblem} =
    LEProblem{ODEProblem}(phase_prob, num_attr; kwargs...)
LEProblem(phase_prob::DEP, num_attr; kwargs...) where {DEP <:DiscreteProblem} =
    LEProblem{DiscreteProblem}(phase_prob, num_attr; kwargs...)

LEProblem{DEP}(phase_prob::DEP;
               num_attr::Int = error("Positional or keyword argument",
                                     " `num_attr` is required."),
               kwargs...) where {DEP} =
    LEProblem{DEP}(phase_prob, num_attr; kwargs...)

struct Relaxer{LEP} <: AbstractRelaxer
    prob::LEP
    integrator
end

get_dim_lyap(integrator) = size(init_tangent_state(integrator))[2]

"""
    LESolver(integrator; <keyword arguments>)

A type representing the main calculation of Lyapunov Exponents (LE).
This struct holds all temporary state required for LE calculation.
"""
mutable struct LESolver{Intr,
                        V <: AbstractVector,
                        M <: AbstractMatrix,
                        } <: AbstractLESolver{Intr}
    integrator::Intr
    num_attr::Int
    series::OnlineStats.Series
    main_stat::OnlineStatsBase.OnlineStat
    inst_exponents::V
    num_orth::Int
    phase_state::V
    tangent_state::M
    Q::M
    R::M
    sign_R::Vector{Bool}
    # TODO: Make sure that they result in concrete types

    function LESolver(
            integrator::Intr, num_attr::Int;
            phase_state::V = init_phase_state(integrator),
            tangent_state::M = init_tangent_state(integrator),
            main_stat = VecMean,
            add_stats = [],
            ) where {Intr,
                     V <: AbstractVector,
                     M <: AbstractMatrix}

        dim_lyap = get_dim_lyap(integrator)
        @assert dim_lyap <= length(phase_state)
        @assert size(tangent_state)[1] == length(phase_state)

        num_orth = 0
        if ! isa(main_stat, OnlineStats.OnlineStat)
            main_stat = main_stat(dim_lyap)
        end
        series = OnlineStats.Series(main_stat, add_stats...)
        inst_exponents = zeros(eltype(phase_state), dim_lyap)
        Q = similar(tangent_state)
        R = similar(tangent_state, (0, 0))  # dummy
        sign_R = Array{Bool}(dim_lyap)
        new{Intr, V, M}(
            integrator,
            num_attr,
            series,
            main_stat,
            inst_exponents,
            num_orth,
            phase_state,
            tangent_state,
            Q,
            R,
            sign_R,
        )
    end
end

"""
    MLESolver(integrator; <keyword arguments>)

A type representing the main calculation of Maximum Lyapunov Exponents
(MLE).  This struct holds all temporary state required for it.
"""
mutable struct MLESolver{Intr,
                         T <: Real,
                         V <: AbstractVector,
                         M <: AbstractMatrix,
                         } <: AbstractLESolver{Intr}
    integrator::Intr
    num_attr::Int
    exponent::T
    inst_exponent::T
    num_orth::Int
    phase_state::V
    tangent_state::M
    # TODO: Make sure that they result in concrete types

    function MLESolver(
            integrator::Intr, num_attr::Int;
            phase_state::V = init_phase_state(integrator),
            tangent_state::M = init_tangent_state(integrator),
            ) where {Intr,
                     V <: AbstractVector,
                     M <: AbstractMatrix}

        if size(tangent_state) != (length(phase_state), 1)
            error("tangent_state must be an array of",
                  " size ($(length(phase_state)), 1);",
                  " given: $(size(tangent_state))")
        end

        num_orth = 0
        T = eltype(phase_state)
        exponent = T(0)
        new{Intr, T, V, M}(
            integrator,
            num_attr,
            exponent,
            exponent,
            num_orth,
            phase_state,
            tangent_state,
        )
    end
end

Base.show(io::IO, solver::LESolver) =
    print(io,
          "#Orth.: ", solver.num_orth, ", ",
          "LEs: ", lyapunov_exponents(solver))
