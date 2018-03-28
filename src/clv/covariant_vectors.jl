module CovariantVectors

# Re-export methods from DifferentialEquations extended here:
export init, solve, step!
import DifferentialEquations: init, solve, step!

using ..Stages: AbstractSource, AbstractStage, finish_if_not!,
    StageIterator, StageState, StagedSolver, goto!
import ..Stages: current_result, stage_index, record!

using ..LyapunovExponents: LEProblem, LESolver, dimension, is_semi_unitary,
    default_Q0, get_tangent_prob, validate_tangent_prob, time_type, ceil_if
import ..LyapunovExponents:
    phase_state, TangentRenormalizer, get_integrator

include("utils.jl")
include("clv_problem.jl")
include("core_stages.jl")
include("clv_solver.jl")
include("accessors.jl")

end

using .CovariantVectors: goto!

using .CovariantVectors:
    CLVProblem, CLVSolver, CLV, backward_dynamics!, forward_dynamics!,
    indexed_backward_dynamics!, indexed_forward_dynamics!
export
    CLVProblem, CLVSolver, CLV, backward_dynamics!, forward_dynamics!,
    indexed_backward_dynamics!, indexed_forward_dynamics!
