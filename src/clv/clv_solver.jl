mutable struct CLVSolution
    # TODO: make those types more general
    G_history::Vector{Matrix{Float64}}
    R_history::Vector{UTM}
    C_history::Vector{UTM}

    CLVSolution() = new()
end

"""
    record!(sol::CLVSolution, stage::AbstractStage)

Record solution of the computation `stage` to solution object `sol`.
"""
function record!(::CLVSolution, ::AbstractStage) end

mutable struct CLVSolver
    prob::CLVProblem
    iter::StageIterator
    state::StageState
    sol::CLVSolution
end

function CLVSolver(prob::CLVProblem)
    stage_types = [
        # PhaseRelaxer,
        ForwardRelaxer,
        ForwardDynamics,
        BackwardRelaxer,
        BackwardDynamics,
    ]
    return CLVSolver(prob, stage_types)
end

function CLVSolver(prob::CLVProblem, stage_types)
    sol = CLVSolution()
    iter = StageIterator(
        prob,
        stage_types,
        (prob,),  # second argument to each `stage_types`
    )
    state = start(iter)
    return CLVSolver(prob, iter, state, sol)
end

CLVSolver(prob::LEProblem) = CLVSolver(CLVProblem(prob))

function advance!(solver::CLVSolver)
    if done(solver.iter, solver.state)
        return nothing
    end
    finish_if_not!(solver.state.stage)
    record!(solver.sol, solver.state.stage)
    _, solver.state = next(solver.iter, solver.state)
    return solver.state.stage
end

function get_last_stage(solver::CLVSolver)
    while advance!(solver) !== nothing end
    return solver.state.stage
end

function solve!(solver::CLVSolver)
    if done(solver.iter, solver.state) && is_finished(solver.state.stage)
        error("No further computation is required.")
        # Should it be just a no-op?
    end
    finish!(get_last_stage(solver))
    return solver
end


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

init(prob::CLVProblem) = CLVSolver(prob)

"""
    forward_dynamics!(solver::CLVSolver) :: ForwardDynamics

Solve the CLV problem up to the forward dynamics stage and return an
iterator to step through the forward dynamics.
"""
forward_dynamics!(solver::CLVSolver) = goto!(solver, ForwardDynamics)

"""
    backward_dynamics!(solver::CLVSolver) :: BackwardDynamics

Solve the CLV problem up to the (final) backward dynamics stage and
return an iterator to step through the backward dynamics.

### Example
```julia
angles = [acos(abs(dot(C[:, 1], C[:, 2]))) * 2 / π for C
          in backward_dynamics!(solver)]
```
"""
backward_dynamics!(solver::CLVSolver) = goto!(solver, BackwardDynamics)
