module Stages

using DiffEqBase: step!

abstract type AbstractStage end
abstract type AbstractComputationStage <: AbstractStage end

"""
    finish!(stage::AbstractStage) :: AbstractStage

Finish whatever the computation `stage` has to do.
"""
function finish! end

function finish!(stage::AbstractComputationStage)
    while ! is_finished(stage)
        step!(stage)
    end
    return stage
end

"""
    is_finished(stage::AbstractStage) :: Bool

Ask if `stage`'s computation is done.
"""
function is_finished end

function finish_if_not!(stage::AbstractStage)
    is_finished(stage) || finish!(stage)
    return stage
end


# Default iterator interface ("mix-in")
current_result(stage::AbstractComputationStage) = nothing

Base.start(stage::AbstractComputationStage) = nothing
Base.done(stage::AbstractComputationStage, _state) = is_finished(stage)
@inline function Base.next(stage::AbstractComputationStage, _state)
    step!(stage)
    return (current_result(stage), nothing)
end
Base.length(stage::AbstractComputationStage) = stage_length(stage)

# TODO: use Base.length directly
function stage_length end


struct StageIterator
    source::AbstractStage
    stage_types::Vector{Type{<: AbstractComputationStage}}
    args::Tuple
end

struct StageState
    i::Int
    stage::AbstractStage
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

end
