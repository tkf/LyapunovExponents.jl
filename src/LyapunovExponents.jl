module LyapunovExponents

using Requires

include("utils.jl")
include("core.jl")
include("coevolve.jl")
include("continuous_exponents.jl")
include("interface.jl")
include("examples.jl")

@require RecipesBase include("plotting.jl")

end # module
