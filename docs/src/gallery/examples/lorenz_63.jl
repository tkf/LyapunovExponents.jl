using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.lorenz_63(num_attr=1000); progress=1)
plot(demo)
