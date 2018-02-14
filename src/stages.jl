module Stages

import DiffEqBase: solve!
using DiffEqBase: step!

abstract type Stageable end
abstract type AbstractSource <: Stageable end
abstract type AbstractStage <: Stageable end

"""
    finish!(stage::Stageable) :: Stageable

Finish whatever the computation `stage` has to do.
"""
function finish! end

function finish!(stage::AbstractStage)
    while ! is_finished(stage)
        step!(stage)
    end
    return stage
end

"""
    is_finished(stage::Stageable) :: Bool

Ask if `stage`'s computation is done.
"""
function is_finished end

function finish_if_not!(stage::Stageable)
    is_finished(stage) || finish!(stage)
    return stage
end

finish!(::AbstractSource) = nothing
is_finished(::AbstractSource) = true


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

function record_finished!(solver::StagedSolver)
    finish_if_not!(solver.state.stage)
end

function advance!(solver::StagedSolver)
    if done(solver.iter, solver.state)
        return nothing
    end
    record_finished!(solver)
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
    record_finished!(solver)
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
