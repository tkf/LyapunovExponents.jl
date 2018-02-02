module LyapunovExponentsTests

using LyapunovExponents
using Base.Test

using Plots
unicodeplots()  # to make tests faster

tic()
include("test_utils.jl")
include("test_online_stats.jl")
include("test_smoke.jl")
include("test_examples.jl")
toc()

end
