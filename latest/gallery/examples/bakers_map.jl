using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.bakers_map(); progress=1)
report(demo)
plt = plot(demo)
