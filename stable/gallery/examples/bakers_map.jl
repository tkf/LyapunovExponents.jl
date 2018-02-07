using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.bakers_map(); progress=1)
plot(demo)
