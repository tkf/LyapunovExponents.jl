using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.van_der_pol(); progress=1)
plot(demo)
