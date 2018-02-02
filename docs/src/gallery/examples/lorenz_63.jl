using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.lorenz_63(); progress=1)
plot(demo)
