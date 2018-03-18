__precompile__()
module LyapunovExponents

export LEProblem, ContinuousLEProblem, DiscreteLEProblem, lyapunov_exponents,
    phase_state, report

# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
import DiffEqBase: init, step!
using DiffEqBase: solve!, set_u!
# solve and solve! are imported/extended in stages.jl

include("documents.jl")
include("stages.jl")
include("utils.jl")
include("online_stats.jl")
include("types.jl")
include("core.jl")
include("terminators.jl")
include("coevolve.jl")
include("continuous_exponents.jl")
include("discrete_exponents.jl")
include("tangent_utils.jl")
include("clv/covariant_vectors.jl")
include("interface.jl")
include("reports.jl")
include("examples/examples.jl")
include("test.jl")
include("diffeq_hack.jl")
include("plotting.jl")

end # module
