using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.henon_map(num_attr=1000); progress=1)
plot(demo)
