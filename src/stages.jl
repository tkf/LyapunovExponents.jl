module Stages

import DiffEqBase: solve!
using DiffEqBase: step!

abstract type Stageable end
abstract type AbstractSource <: Stageable end
abstract type AbstractStage <: Stageable end

"""
    is_finished(stage::Stageable) :: Bool

Ask if `stage`'s computation is done.
"""
function is_finished end

"""
    finish!(stage::Stageable) :: Stageable

Finish whatever the computation `stage` has to do.
"""
function finish!(stage::AbstractStage)
    while ! is_finished(stage)
        step!(stage)
    end
    return stage
end

function finish_if_not!(stage::Stageable)
    is_finished(stage) || finish!(stage)
    return stage
end

finish!(::AbstractSource) = nothing
is_finished(::AbstractSource) = true

# Some "mix-in" methods for generating default `is_finished`:
stage_index(stage::AbstractStage) = stage.i
is_finished(stage::AbstractStage) = stage_index(stage) >= length(stage)
# TODO: maybe this shouldn't be defined for generic stage.  Maybe
# define AbstractIndexedStage?


# Utility method. Define it and call it from step! for type-based
# recording injection.
"""
    record!(stage, Val{key})

Record some info (labeled by `key`) into `stage`.
"""
record!(::AbstractStage, ::Any) = nothing
# See: [[./clv/core_stages.jl::function record!]]


# Default iterator interface ("mix-in")
current_result(stage::AbstractStage) = nothing

Base.start(stage::AbstractStage) = nothing
Base.done(stage::AbstractStage, _state) = is_finished(stage)
@inline function Base.next(stage::AbstractStage, _state)
    step!(stage)
    return (current_result(stage), nothing)
end


struct StageIterator
    source::Stageable
    stage_types::Vector{Type{<: AbstractStage}}
    args::Tuple
end

struct StageState
    i::Int
    stage::Stageable
end


Base.start(iter::StageIterator) = StageState(1, iter.source)

function Base.next(iter::StageIterator, state::StageState)
    next_state = StageState(
        state.i + 1,
        iter.stage_types[state.i](
            finish_if_not!(state.stage),
            iter.args...))
    return (next_state.stage, next_state)
end

Base.done(iter::StageIterator, state::StageState) =
    length(iter.stage_types) < state.i


is_reachable(iter::StageIterator, stage_type::Type, i::Int = 1) =
    any(t <: stage_type for t in iter.stage_types[i:end])

is_reachable(iter::StageIterator, stage_type::Type, state::StageState) =
    is_reachable(iter, stage_type, state.i)


"""
Stateful solver based on [`StageIterator`](@ref).
"""
mutable struct StagedSolver{P, S}
    prob::P
    iter::StageIterator
    state::StageState
    sol::S

    function StagedSolver(prob::P, sol::S,
                          stage_types::AbstractVector, args::Tuple,
                          ) where {P, S}
        iter = StageIterator(
            prob,
            stage_types,
            args,  # additional arguments to each `stage_types`
        )
        state = start(iter)
        return new{P, S}(prob, iter, state, sol)
    end
end

function advance!(solver::StagedSolver)
    if done(solver.iter, solver.state)
        return nothing
    end
    finish_if_not!(solver.state.stage)
    _, solver.state = next(solver.iter, solver.state)
    return solver.state.stage
end

function get_last_stage!(solver::StagedSolver)
    while advance!(solver) !== nothing end
    return solver.state.stage
end

function solve!(solver::StagedSolver)
    if done(solver.iter, solver.state) && is_finished(solver.state.stage)
        error("No further computation is required.")
        # Should it be just a no-op?
    end
    get_last_stage!(solver)
    finish_if_not!(solver.state.stage)
    return solver
end


"""
    goto!(solver::StagedSolver, stage_type::Type{T}) :: T

Advance the `solver` up to the stage of type `stage_type` and return
it.
"""
function goto!(solver::StagedSolver, stage_type::Type)
    @assert is_reachable(solver.iter, stage_type, solver.state)
    if solver.state.stage isa stage_type
        return solver.state.stage
    end
    while true
        stage = advance!(solver)
        if stage isa stage_type
            return stage
        end
    end
    error("Unreachable!")
end

end
