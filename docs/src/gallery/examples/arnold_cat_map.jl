using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.arnold_cat_map(); progress=1)
report(demo)
plt = plot(demo)
