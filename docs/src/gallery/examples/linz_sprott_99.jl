using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.linz_sprott_99(num_attr=1000); progress=1)
plot(demo)
