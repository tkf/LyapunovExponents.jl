using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.linz_sprott_99(); progress=1)
plt = plot(demo)
