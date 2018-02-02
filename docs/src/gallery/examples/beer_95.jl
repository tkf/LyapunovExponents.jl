using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.beer_95(num_attr=1000); progress=1)
plot(demo)
