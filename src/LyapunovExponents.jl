module LyapunovExponents

export LEProblem, ContinuousLEProblem, DiscreteLEProblem, lyapunov_exponents,
    phase_state

using Requires

include("stages.jl")
include("utils.jl")
include("online_stats.jl")
include("types.jl")
include("core.jl")
include("coevolve.jl")
include("continuous_exponents.jl")
include("discrete_exponents.jl")
include("recording.jl")
include("tangent_utils.jl")
include("clv/covariant_vectors.jl")
include("interface.jl")
include("examples/examples.jl")
include("test.jl")
include("diffeq_hack.jl")

@require RecipesBase include("plotting.jl")

end # module
