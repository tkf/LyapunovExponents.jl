module Stages

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

end
