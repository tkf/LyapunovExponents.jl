module LyapunovExponents

export ContinuousLEProblem, DiscreteLEProblem, lyapunov_exponents

using Requires

include("utils.jl")
include("types.jl")
include("core.jl")
include("coevolve.jl")
include("continuous_exponents.jl")
include("discrete_exponents.jl")
include("recording.jl")
include("interface.jl")
include("examples.jl")

@require RecipesBase include("plotting.jl")

end # module
