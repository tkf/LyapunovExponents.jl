using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.henon_map(); progress=1)
plt = plot(demo)
