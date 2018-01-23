module LyapunovExponentsTests

using LyapunovExponents
using Base.Test

tic()
include("test_examples.jl")
include("test_ode.jl")
toc()

end
