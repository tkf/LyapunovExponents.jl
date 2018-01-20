function lyapunov_exponents(phase_dynamics!,
                            u0,
                            tspan,
                            num_attr::Integer;
                            progress=-1,
                            kwargs...)
    prob = ContinuousLEProblem(phase_dynamics!, u0, tspan, kwargs...)
    solver = init(prob; progress=progress)
    solve!(solver, num_attr; progress=progress)
    lyapunov_exponents(solver)
end
