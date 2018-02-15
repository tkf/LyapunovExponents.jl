module LyapunovExponents

export LEProblem, ContinuousLEProblem, DiscreteLEProblem, lyapunov_exponents,
    phase_state

# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
import DiffEqBase: init, step!
using DiffEqBase: solve, solve!
# solve and solve! are imported/extended in stages.jl

using Requires

include("stages.jl")
include("utils.jl")
include("online_stats.jl")
include("types.jl")
include("core.jl")
include("coevolve.jl")
include("continuous_exponents.jl")
include("discrete_exponents.jl")
include("tangent_utils.jl")
include("clv/covariant_vectors.jl")
include("interface.jl")
include("examples/examples.jl")
include("test.jl")
include("diffeq_hack.jl")

@require RecipesBase include("plotting.jl")

end # module
