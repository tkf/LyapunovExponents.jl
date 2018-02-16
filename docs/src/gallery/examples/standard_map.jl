using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.standard_map(); progress=1)
plt = plot(demo)
