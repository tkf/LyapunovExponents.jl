using LyapunovExponents

function get_angles(solver)
    return [acos(abs(dot(C[:, 1], C[:, 2]))) * 2 / π for C
            in goto!(solver, BackwardDynamics)]
end

henon_demo = LyapunovExponents.henon_map(num_attr=100000)
henon_demo.prob :: LEProblem
henon_prob = CLVProblem(henon_demo.prob)  # convert it to CLVProblem
henon_angles = @time get_angles(init(henon_prob))

lozi_prob = CLVProblem(LyapunovExponents.lozi_map(num_attr=100000).prob)
lozi_angles = @time get_angles(init(lozi_prob))

using Plots
plt = plot(xlabel="Angle [pi/2 rad]", ylabel="Density", legend=:topleft)
stephist!(plt, henon_angles,
          bins=1000, normalize=true, linecolor=1, label="Henon")
stephist!(twinx(plt), lozi_angles,
          bins=1000, normalize=true, linecolor=2, label="Lozi")
plt
