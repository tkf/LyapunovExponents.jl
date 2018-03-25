module LyapunovExponentsTests

print_with_color(:blue, let
    message = " Warnings below (if any) are fine. "
    margin = (displaysize(STDOUT)[2] - length(message)) ÷ 2
    ("=" ^ margin) * message * ("=" ^ margin)
end)
println()
flush(STDOUT)
import Plots
import ForwardDiff
import DifferentialEquations
import OnlineStats
flush(STDOUT)
flush(STDERR)
print_with_color(:blue, "=" ^ displaysize(STDOUT)[2])
println()

using LyapunovExponents
using Base.Test

tic()
include("test_utils.jl")
include("test_online_stats.jl")
include("test_smoke.jl")
include("test_ui.jl")
include("test_examples.jl")
include("test_terminators.jl")
include("test_clv.jl")
include("test_null_clv.jl")
toc()

end
