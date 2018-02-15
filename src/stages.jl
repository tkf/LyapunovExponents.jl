module Stages

import DiffEqBase: solve!, solve
using DiffEqBase: step!, init

using ProgressMeter: ProgressMeter, Progress

abstract type Stageable end
abstract type AbstractSource <: Stageable end
abstract type AbstractStage <: Stageable end

abstract type ExecPolicy end
struct NullPolicy <: ExecPolicy end
struct WithProgress{P} <: ExecPolicy
    progress_factory::P
    WithProgress(progress_factory::P) where P = new{P}(progress_factory)
end

get_exec_policy(progress::Real) =
    if progress > 0
        WithProgress((n, args...) -> Progress(n, progress, args...))
    else
        NullPolicy()
    end

pre_loop(::AbstractStage, ::NullPolicy) = nothing
post_loop(::AbstractStage, ::NullPolicy, ::Any) = nothing
post_step(::AbstractStage, ::NullPolicy, ::Any) = nothing
function remaining_steps end

"""
    is_finished(stage::Stageable) :: Bool

Ask if `stage`'s computation is done.
"""
function is_finished end

"""
    finish!(stage::Stageable) :: Stageable

Finish whatever the computation `stage` has to do.
"""
function finish!(stage::AbstractStage, policy::ExecPolicy = NullPolicy())
    policy_state = pre_loop(stage, policy)
    while ! is_finished(stage)
        step!(stage)
        policy_state = post_step(stage, policy, policy_state)
    end
    post_loop(stage, policy, policy_state)
    return stage
end

shortname(::T) where T = string(T.name.name)

pre_loop(stage::AbstractStage, policy::WithProgress) =
    policy.progress_factory(remaining_steps(stage), shortname(stage))

function post_step(::AbstractStage, ::WithProgress, prog::Progress)
    ProgressMeter.next!(prog)
    return prog
end

post_loop(::AbstractStage, ::WithProgress, prog::Progress) =
    ProgressMeter.finish!(prog)

function finish_if_not!(stage::Stageable, policy::ExecPolicy = NullPolicy())
    is_finished(stage) || finish!(stage, policy)
    return stage
end

finish!(::AbstractSource, ::ExecPolicy) = nothing
is_finished(::AbstractSource) = true

# Some "mix-in" methods for generating default `is_finished`:
stage_index(stage::AbstractStage) = stage.i
is_finished(stage::AbstractStage) = stage_index(stage) >= length(stage)
remaining_steps(stage::AbstractStage) = length(stage) - stage_index(stage)
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

function advance!(solver::StagedSolver, policy::ExecPolicy = NullPolicy())
    if done(solver.iter, solver.state)
        return nothing
    end
    finish_if_not!(solver.state.stage, policy)
    _, solver.state = next(solver.iter, solver.state)
    return solver.state.stage
end

function get_last_stage!(solver::StagedSolver,
                         policy::ExecPolicy = NullPolicy())
    while advance!(solver, policy) !== nothing end
    return solver.state.stage
end

function solve!(solver::StagedSolver; progress=-1)
    if done(solver.iter, solver.state) && is_finished(solver.state.stage)
        error("No further computation is required.")
        # Should it be just a no-op?
    end
    policy = get_exec_policy(progress)
    get_last_stage!(solver, policy)
    finish_if_not!(solver.state.stage, policy)
    return solver
end

solve(prob::AbstractSource; progress=-1, kwargs...) =
    solve!(init(prob; kwargs...); progress=progress).sol
# This requires `init(prob)` to be defined elsewhere and to return a
# StagedSolver.  Not sure defining `solve` by default is a good idea.
# But repeating this definition all over the places is also not good.

"""
    goto!(solver::StagedSolver, stage_type::Type{T}) :: T

Advance the `solver` up to the stage of type `stage_type` and return
it.
"""
function goto!(solver::StagedSolver, stage_type::Type; progress=-1)
    @assert is_reachable(solver.iter, stage_type, solver.state)
    if solver.state.stage isa stage_type
        return solver.state.stage
    end
    policy = get_exec_policy(progress)
    while true
        stage = advance!(solver, policy)
        if stage isa stage_type
            return stage
        end
    end
    error("Unreachable!")
end

end
