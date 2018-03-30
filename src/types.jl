using DiffEqBase: DEProblem, ODEProblem, DiscreteProblem, remake,
    parameterless_type
using .Stages: Stageable, AbstractSource, AbstractStage,
    StageIterator, StageState, StagedSolver, goto!
import .Stages: record!, current_result, is_finished

function time_type(is_discrete; times...)
    TT = promote_type((typeof(v) for (_, v) in times)...)

    if is_discrete && ! (TT <: Integer)
        msg = "Non integer time(s):"
        for (name, val) in times
            if ! (val isa Integer)
                msg = "$msg\n  $name::$(typeof(val)) = $val"
            end
        end
        error(msg)
    end

    return TT
end

"""
    LEProblem(phase_prob, t_attr; <keyword arguments>)

# Arguments
- `phase_prob`: Phase space dynamics represented in the form of
  `ODEProblem` or `DiscreteProblem` from DifferentialEquations.jl.
  `phase_prob.tspan` is ignored.
- `t_attr::Real`: Simulated time on the (presumed) attractor
  used for the solver.  It is used for computation of Lyapunov
  Exponents.  Roughly `t_attr / t_renorm` instantaneous exponents are
  sampled.
- `t_tran::Real`: Number of iterations to throw away to get rid of the
  effect from the transient dynamics.
- `t_renorm::Real`: Interval between orthonormalizations (renormalizations).
- `dim_lyap::Integer`: Number of Lyapunov exponents to be calculated.
  Default to the full system dimension.
- `Q0::Array`: The initial guess of the Gram-Schmidt "Lyapunov vectors".
- `tangent_dynamics::Function`: A vector field for solving phase space
   evolution *and* tangent space evolution together.  If this is not
   provided, `tangent_dynamics` is derived from `phase_prob.f`.  See
   also [`PhaseTangentDynamics`](@ref).
"""
struct LEProblem{PPr, TPr, TT} <: AbstractSource
    phase_prob::PPr
    tangent_prob::TPr
    t_tran::TT
    t_attr::TT
    t_renorm::TT

    function LEProblem{PPr, TPr}(
            phase_prob::PPr,
            tangent_prob::TPr,
            t_attr::Real;
            t_renorm::Real = 1,
            t_tran::Real = 1,
    ) where {PPr, TPr}

        TT = time_type((PPr <: DiscreteProblem);
                       t_tran = t_tran,
                       t_attr = t_attr,
                       t_renorm = t_renorm)

        new{PPr, TPr, TT}(
            phase_prob, tangent_prob,
            t_tran, t_attr, t_renorm)
    end
end

make_tangent_prob(phase_prob, tangent_dynamics, u0) =
    remake(phase_prob,
           f = tangent_dynamics,
           u0 = u0)

function validate_tangent_prob(phase_prob::DEProblem,
                               tangent_prob,
                               tangent_dynamics,
                               Q0)
    if tangent_prob !== nothing
        if tangent_dynamics !== nothing
            error("Both tangent_prob and tangent_dynamics are given. ",
                  "Only at most one of them can be specified.")
        end
    else
        x0 = phase_prob.u0
        u0 = similar(x0, (length(x0), size(Q0, 2) + 1))
        u0[:, 1] = x0  # ignored anyway
        u0[:, 2:end] = Q0
        if tangent_dynamics == nothing
            tangent_dynamics = fdiff_tangent_dynamics(phase_prob.f, u0)
        end
        tangent_prob = make_tangent_prob(phase_prob, tangent_dynamics, u0)
    end
    if ! is_semi_unitary(@view tangent_prob.u0[:, 2:end])
        error("Columns in Q0 are not orthonormal.\n",
              "Q0 = tangent_prob.u0[:, 2:end]")
    end
    return tangent_prob
end

function LEProblem(phase_prob::PPr;
                   t_attr::Real = error("Keyword argument",
                                        " `t_attr` is required."),
                   tangent_dynamics = nothing,
                   tangent_prob::Union{Void, DEProblem} = nothing,
                   dim_lyap = dimension(phase_prob),
                   Q0 = default_Q0(phase_prob, dimension(phase_prob), dim_lyap),
                   kwargs...) where {PPr <: DEProblem}
    tangent_prob = validate_tangent_prob(
        phase_prob, tangent_prob, tangent_dynamics, Q0)
    LEProblem{parameterless_type(PPr),
              parameterless_type(typeof(tangent_prob))}(
        phase_prob, tangent_prob, t_attr;
        kwargs...)
end

mutable struct PhaseRelaxer{Intr, T} <: AbstractStage
    integrator::Intr
    t_tran::T
    i::Int

    PhaseRelaxer(integrator::Intr, t_tran::T) where {Intr, T} =
        new{Intr, T}(integrator, t_tran, 0)
end


mutable struct LESolution{RecOS, RecFTLE, SS, MS, VV, CH}
    series::SS
    main_stat::MS
    ftle_history::VV
    convergence::CH
    num_orth::Int
    # TODO: Bundle ftle_history and num_orth in a struct.
    converged::Bool
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
    convergence = if RecOS && RecFTLE
        ConvergenceHistory(dim_lyap)
    end
    args = (series, main_stat, ftle_history, convergence)
    return LESolution{RecOS, RecFTLE, map(typeof, args)...}(
        args...,
        0,                      # num_orth
        false,                  # converged
    )
end

NullLESolution() = LESolution(0; main_stat = nothing)

# TODO: define record!(sol::LESolution, ftle)

"""
    get_dim_lyap(thing)::Int
"""
function get_dim_lyap end

abstract type Terminator end

abstract type AbstractRenormalizer{S} <: AbstractStage end

"""
    TangentRenormalizer(integrator; <keyword arguments>)

A type representing the main calculation of Lyapunov Exponents (LE).
This struct holds all temporary state required for LE calculation.
"""
mutable struct TangentRenormalizer{S <: LESolution,
                                   TMNR <: Terminator,
                                   Intr, T, V, M,
                                   } <: AbstractRenormalizer{S}
    integrator::Intr
    t_attr::T
    t_renorm::T
    inst_exponents::V
    i::Int
    phase_state::V
    tangent_state::M
    Q::M
    R::M
    sign_R::Vector{Bool}
    # TODO: Make sure that they result in concrete types

    sol::S
    tmnr::TMNR
    is_finished::Bool

    function TangentRenormalizer(
            integrator::Intr, t_attr, t_renorm,
            sol::S = NullLESolution(),
            tmnr::TMNR = NullTerminator();
            phase_state::V = init_phase_state(integrator),
            tangent_state::M = init_tangent_state(integrator),
            ) where {S, Intr,
                     TMNR <: Terminator,
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

        T = promote_type(map(typeof, (t_attr, t_renorm))...)
        return new{S, TMNR, Intr, T, V, M}(
            integrator,
            t_attr,
            t_renorm,
            inst_exponents,
            0,  # i
            phase_state,
            tangent_state,
            Q,
            R,
            sign_R,
            sol,
            tmnr,
            false,  # is_finished
        )
    end
end

mutable struct MLERenormalizer{S <: LESolution,
                               TMNR <: Terminator,
                               Intr,
                               T <: Real,
                               E <: Real,
                               V <: AbstractVector,
                               M <: AbstractMatrix,
                               } <: AbstractRenormalizer{S}
    integrator::Intr
    t_attr::T
    t_renorm::T
    exponent::E
    inst_exponent::E
    i::Int
    phase_state::V
    tangent_state::M
    # TODO: Make sure that they result in concrete types

    sol::S
    tmnr::TMNR
    is_finished::Bool

    function MLERenormalizer(
            integrator::Intr, t_attr, t_renorm,
            sol::S = NullLESolution(),
            tmnr::TMNR = NullTerminator();
            phase_state::V = init_phase_state(integrator),
            tangent_state::M = init_tangent_state(integrator),
            ) where {S, Intr,
                     TMNR <: Terminator,
                     V <: AbstractVector,
                     M <: AbstractMatrix}

        if size(tangent_state) != (length(phase_state), 1)
            error("tangent_state must be an array of",
                  " size ($(length(phase_state)), 1);",
                  " given: $(size(tangent_state))")
        end

        E = eltype(phase_state)
        exponent = E(0)
        T = promote_type(map(typeof, (t_attr, t_renorm))...)
        new{S, TMNR, Intr, T, E, V, M}(
            integrator,
            t_attr,
            t_renorm,
            exponent,
            exponent,
            0,  # i
            phase_state,
            tangent_state,
            sol,
            tmnr,
            false,  # is_finished
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
                  renormalizer = if get_dim_lyap(prob) == 1
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
                  integrator_options = [],
                  terminator = false,
                  terminator_options = [],
                  kwargs...)
    sol = LESolution(get_dim_lyap(prob);
                     num_attr = if record
                         ceil(Int, prob.t_attr / prob.t_renorm)
                     else
                         nothing
                     end,
                     kwargs...)
    tmnr = Terminator(terminator, prob, sol; terminator_options...)
    # Additional arguments to each `stage_types`:
    args = (prob, sol, tmnr, integrator_options)
    return StagedSolver(prob, sol, stage_types, args)
end

init(prob::LEProblem; kwargs...) = LESolver(prob; kwargs...)

num_orth(stage::AbstractRenormalizer) = stage.i
num_orth(stage::Stageable) = 0
num_orth(solver::LESolver) = num_orth(solver.state.stage)

Base.show(io::IO, solver::LESolver) =
    print(io,
          "#Orth.: ", num_orth(solver), ", ",
          "LEs: ", lyapunov_exponents(solver))
