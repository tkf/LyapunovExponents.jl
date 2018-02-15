using DiffEqBase: ODEProblem, DiscreteProblem
using .Stages: AbstractSource, AbstractStage,
    StageIterator, StageState, StagedSolver, goto!
import .Stages: record!

"""
    LEProblem(phase_prob, num_attr; <keyword arguments>)
    LEProblem(phase_prob; num_attr, <keyword arguments>)

# Arguments
- `phase_prob`: Phase space dynamics represented in the form of
  `ODEProblem` or `DiscreteProblem` from DifferentialEquations.jl.
  `phase_prob.tspan` represents the inter-orthonormalization-interval.
- `num_attr::Integer`: Number of orthonormalizations (or some kind of
  "renormalisation") for calculating Lyapunov Exponents.  In general,
  this is the number of points (considered to be) on the attractor
  used for the solver.  The simulated time of the system for this
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
struct LEProblem{DEP} <: AbstractSource
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

mutable struct PhaseRelaxer{Intr} <: AbstractStage
    integrator::Intr
    num_tran::Int
    i::Int

    PhaseRelaxer(integrator::Intr, num_tran::Int, i::Int = 1) where {Intr} =
        new{Intr}(integrator, num_tran, i)
end

struct LESolution{RecOS, RecFTLE, SS, MS, VV}
    series::SS
    main_stat::MS
    ftle_history::VV
end

const LESolRecOS = LESolution{true}
const LESolRecFTLE = LESolution{RecOS, true} where {RecOS}

function LESolution(dim_lyap::Int;
                    main_stat = VecMean,
                    add_stats = [],
                    history_type = Vector{Float64},
                    num_attr = nothing)
    RecOS = main_stat !== nothing || ! isempty(add_stats)
    RecFTLE = num_attr !== nothing
    if ! isa(main_stat, Union{OnlineStats.OnlineStat, Void})
        main_stat = main_stat(dim_lyap)
    end
    series = if RecOS
        OnlineStats.Series(main_stat, add_stats...)
    end
    ftle_history = if RecFTLE
        [history_type(dim_lyap) for _ in 1:num_attr]
        # TODO: generalize outer array type?
    end
    args = (series, main_stat, ftle_history)
    return LESolution{RecOS, RecFTLE, map(typeof, args)...}(args...)
end

NullLESolution() = LESolution(0; main_stat = nothing)

get_dim_lyap(integrator) = size(init_tangent_state(integrator))[2]

abstract type AbstractRenormalizer{S} <: AbstractStage end

"""
    TangentRenormalizer(integrator; <keyword arguments>)

A type representing the main calculation of Lyapunov Exponents (LE).
This struct holds all temporary state required for LE calculation.
"""
mutable struct TangentRenormalizer{S <: LESolution,
                                   Intr, V, M,
                                   } <: AbstractRenormalizer{S}
    integrator::Intr
    num_attr::Int
    inst_exponents::V
    i::Int
    phase_state::V
    tangent_state::M
    Q::M
    R::M
    sign_R::Vector{Bool}
    # TODO: Make sure that they result in concrete types

    sol::S

    function TangentRenormalizer(
            integrator::Intr, num_attr::Int,
            sol::S = NullLESolution();
            phase_state::V = init_phase_state(integrator),
            tangent_state::M = init_tangent_state(integrator),
            ) where {S, Intr,
                     V <: AbstractVector,
                     M <: AbstractMatrix}

        dim_lyap = get_dim_lyap(integrator)
        dim_phase = length(init_phase_state(integrator))
        @assert dim_lyap <= dim_phase
        @assert size(phase_state) == (dim_phase,)
        @assert size(tangent_state) == (dim_phase, dim_lyap)

        inst_exponents = zeros(eltype(phase_state), dim_lyap)
        Q = similar(tangent_state)
        R = similar(tangent_state, (0, 0))  # dummy
        sign_R = Array{Bool}(dim_lyap)

        return new{S, Intr, V, M}(
            integrator,
            num_attr,
            inst_exponents,
            0,  # i
            phase_state,
            tangent_state,
            Q,
            R,
            sign_R,
            sol,
        )
    end
end

mutable struct MLERenormalizer{S <: LESolution,
                               Intr,
                               T <: Real,
                               V <: AbstractVector,
                               M <: AbstractMatrix,
                               } <: AbstractRenormalizer{S}
    integrator::Intr
    num_attr::Int
    exponent::T
    inst_exponent::T
    i::Int
    phase_state::V
    tangent_state::M
    # TODO: Make sure that they result in concrete types

    sol::S

    function MLERenormalizer(
            integrator::Intr, sol::S, num_attr::Int;
            phase_state::V = init_phase_state(integrator),
            tangent_state::M = init_tangent_state(integrator),
            ) where {S, Intr,
                     V <: AbstractVector,
                     M <: AbstractMatrix}

        if size(tangent_state) != (length(phase_state), 1)
            error("tangent_state must be an array of",
                  " size ($(length(phase_state)), 1);",
                  " given: $(size(tangent_state))")
        end

        T = eltype(phase_state)
        exponent = T(0)
        new{S, Intr, T, V, M}(
            integrator,
            num_attr,
            exponent,
            exponent,
            0,  # i
            phase_state,
            tangent_state,
            sol,
        )
    end
end


"""
    phase_state(stage) :: Vector

Get current phase-space state stored in `stage`.
"""
phase_state(stage::TangentRenormalizer) = stage.phase_state
phase_state(stage::MLERenormalizer) = stage.phase_state


const LESolver = StagedSolver{<: LEProblem, <: LESolution}
const LESolverRecOS = StagedSolver{<: LEProblem, <: LESolRecOS}
const LESolverRecFTLE = StagedSolver{<: LEProblem, <: LESolRecFTLE}

"""
    LESolver(prob::LEProblem; record::Bool = false)

Create a solver object for a [`LEProblem`](@ref).  Record all
finite-time (instantaneous) Lyapunov exponents when `record = true` is
passed.
"""
function LESolver(prob::LEProblem;
                  phase_relaxer = PhaseRelaxer,
                  renormalizer = if prob.dim_lyap == 1
                      MLERenormalizer
                  else
                      TangentRenormalizer
                  end,
                  kwargs...)
    stage_types = [
        phase_relaxer,
        renormalizer,
    ]
    return LESolver(prob, stage_types; kwargs...)
end

function LESolver(prob::LEProblem, stage_types::AbstractVector;
                  record::Bool = false,
                  kwargs...)
    sol = LESolution(prob.dim_lyap;
                     num_attr = record ? prob.num_attr : nothing,
                     kwargs...)
    args = (prob, sol)  # additional arguments to each `stage_types`
    return StagedSolver(prob, sol, stage_types, args)
end

init(prob::LEProblem; kwargs...) = LESolver(prob; kwargs...)


Base.show(io::IO, solver::LESolver) =
    print(io,
          "#Orth.: ", solver.i, ", ",
          "LEs: ", lyapunov_exponents(solver))
