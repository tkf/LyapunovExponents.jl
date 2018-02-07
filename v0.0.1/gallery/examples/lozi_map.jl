using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.lozi_map(); progress=1)
plot(demo)
