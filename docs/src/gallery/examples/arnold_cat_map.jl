using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.arnold_cat_map(); progress=1)
plt = plot(demo)
