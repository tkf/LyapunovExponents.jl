using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.van_der_pol(num_attr=1000); progress=1)
plot(demo)
