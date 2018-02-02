using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.bakers_map(num_attr=1000); progress=1)
plot(demo)
