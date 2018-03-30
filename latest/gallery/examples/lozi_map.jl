using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.lozi_map(); progress=1)
report(demo)
plt = plot(demo)
