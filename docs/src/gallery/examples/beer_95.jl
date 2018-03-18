using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.beer_95(); progress=1)
report(demo)
plt = plot(demo)
