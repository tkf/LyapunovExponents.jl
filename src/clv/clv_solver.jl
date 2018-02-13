using Base: Callable

mutable struct CLVSolver
    prob::CLVProblem
    iter::StageIterator
    state::StageState
    sol::CLVSolution
end

"""
    CLVSolver(prob::CLVProblem; <keyword arguments>)
    CLVSolver(prob::LEProblem; <keyword arguments>)

The preferred and equivalent method to get a solver for a `CLVProblem`
is `init(prob::CLVProblem)`.  Note that `CLVSolver(prob::LEProblem)`
is equivalent to `init(CLVProblem(prob))`.
"""
function CLVSolver(prob::CLVProblem;
                   forward_relaxer::Type = ForwardRelaxer,
                   forward_dynamics::Type = ForwardDynamics,
                   backward_relaxer::Type = BackwardRelaxer,
                   backward_dynamics::Type = BackwardDynamics,
                   record::Vector{Symbol} = Symbol[])
    stage_types = [
        # PhaseRelaxer???
        forward_relaxer,
        forward_dynamics,
        backward_relaxer,
        backward_dynamics,
    ]
    return CLVSolver(prob, stage_types, record)
end

function CLVSolver(prob::CLVProblem, stage_types::AbstractVector,
                   record::Vector{Symbol})
    sol = CLVSolution(prob, record)
    iter = StageIterator(
        prob,
        stage_types,
        (prob, sol),  # additional arguments to each `stage_types`
    )
    state = start(iter)
    return CLVSolver(prob, iter, state, sol)
end

CLVSolver(prob::LEProblem; kwargs...) = CLVSolver(CLVProblem(prob); kwargs...)

function record_finished!(solver::CLVSolver)
    finish_if_not!(solver.state.stage)
end

function advance!(solver::CLVSolver)
    if done(solver.iter, solver.state)
        return nothing
    end
    record_finished!(solver)
    _, solver.state = next(solver.iter, solver.state)
    return solver.state.stage
end

function get_last_stage!(solver::CLVSolver)
    while advance!(solver) !== nothing end
    return solver.state.stage
end

function solve!(solver::CLVSolver)
    if done(solver.iter, solver.state) && is_finished(solver.state.stage)
        error("No further computation is required.")
        # Should it be just a no-op?
    end
    get_last_stage!(solver)
    record_finished!(solver)
    return solver
end


"""
    goto!(solver::CLVSolver, stage_type::Type{T}) :: T

Solve the CLV problem up to the stage of type `stage_type` and return
it.  Instead of directly calling `goto!(solver, ForwardDynamics)` and
`goto!(solver, BackwardDynamics)`, use the convenience methods
[`forward_dynamics!`](@ref) and [`backward_dynamics!`](@ref).
"""
function goto!(solver::CLVSolver, stage_type::Type)
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

init(prob::CLVProblem; kwargs...) = CLVSolver(prob; kwargs...)

"""
    forward_dynamics!(solver::CLVSolver) :: ForwardDynamics

Solve the CLV problem up to the forward dynamics stage and return an
iterator to step through the forward dynamics.
See also: [`backward_dynamics!`](@ref), [`goto!`](@ref).
"""
forward_dynamics!(solver::CLVSolver) = goto!(solver, ForwardDynamics)

"""
    indexed_forward_dynamics!(solver::CLVSolver)
    indexed_forward_dynamics!(stage::ForwardDynamics)

Just a short-hand for `enumerate(forward_dynamics!(solver))`.
It's for symmetry with [`indexed_backward_dynamics!`](@ref).
"""
indexed_forward_dynamics!(solver::CLVSolver) =
    indexed_forward_dynamics!(forward_dynamics!(solver))
indexed_forward_dynamics!(fitr::ForwardDynamics) = enumerate(fitr)

"""
    backward_dynamics!(solver::CLVSolver) :: BackwardDynamics

Solve the CLV problem up to the (final) backward dynamics stage and
return an iterator to step through the backward dynamics.
See also [`forward_dynamics!`](@ref), [`goto!`](@ref).

Note that finishing iteration of the returned iterator does not
finalize all the solver stages (namely, recording to `solver.sol`, if
non-default `backward_dynamics` is used).  In this case,
`solve!(solver)` has to be called after the iteration.

### Example
```julia
angles = [acos(abs(dot(C[:, 1], C[:, 2]))) * 2 / π for C
          in backward_dynamics!(solver)]
```
"""
backward_dynamics!(solver::CLVSolver) = goto!(solver, BackwardDynamics)

"""
    indexed_backward_dynamics!(solver::CLVSolver)
    indexed_backward_dynamics!(stage::BackwardDynamics)

It is equivalent to `zip(some_counter, backward_dynamics!(solver))`
where `some_counter` is an iterator over integers such that the same
indices returned by [`indexed_forward_dynamics!`](@ref) indicate that
those events are at the same time point of the backward and forward
passes.  This is useful when combining matrices [`CLV.G`](@ref) and
[`CLV.C`](@ref) to obtain the CLV in the (original) tangent space.

For such example, see:
[Covariant Lyapunov vectors on the Lorenz attractor](@ref)
in the online manual.

Note that `indexed_backward_dynamics!(solver::CLVSolver)` is
equivalent to
```julia
backward = backward_dynamics!(solver)
indexed_backward_dynamics!(backward)
```
Separately calling [`backward_dynamics!`](@ref) is useful when the
quantities other than [`CLV.G`](@ref) (e.g., [`CLV.R`](@ref)) are
required.

See also [`indexed_forward_dynamics!`](@ref).
"""
indexed_backward_dynamics!(solver::CLVSolver) =
    indexed_backward_dynamics!(backward_dynamics!(solver))

function indexed_backward_dynamics!(bitr::BackwardDynamics;
                                    until::Int = 1)
    return zip(length(bitr)-1:-1:until, bitr)
end
