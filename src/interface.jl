"""
    lyapunov_exponents(phase_dynamics!, u0, tspan; <keyword arguments>)

Calculate Lyapunov exponents of a dynamical system.
"""
function lyapunov_exponents(phase_dynamics!,
                            u0,
                            tspan,
                            num_attr::Integer;
                            discrete=false,
                            progress=-1,
                            kwargs...)
    if discrete
        prob = ContinuousLEProblem(phase_dynamics!, u0, tspan, kwargs...)
    else
        prob = DiscreteLEProblem(phase_dynamics!, u0, tspan, kwargs...)
    end
    solver = init(prob; progress=progress)
    solve!(solver, num_attr; progress=progress)
    lyapunov_exponents(solver)
end
