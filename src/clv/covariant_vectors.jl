const get_le_solver = get_solver

module CovariantVectors

# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
import DifferentialEquations: init, solve, solve!, step!

using ..LyapunovExponents: LEProblem, LESolver, dimension, is_semi_unitary,
    default_Q0, get_tangent_dynamics
import ..LyapunovExponents: phase_tangent_state, get_tangent_prob,
    get_le_solver, phase_state

include("utils.jl")
include("stage_base.jl")
include("clv_problem.jl")
include("core_stages.jl")
include("clv_solver.jl")
include("recording_stages.jl")
include("accessors.jl")

end

using .CovariantVectors: goto!

using .CovariantVectors:
    CLVProblem, CLVSolver, CLV, backward_dynamics!, forward_dynamics!
export
    CLVProblem, CLVSolver, CLV, backward_dynamics!, forward_dynamics!
