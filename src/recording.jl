mutable struct LERecordingSolver{Intr} <: AbstractLESolver{Intr}
    solver
    exponents_history
    i_hist
end

function LERecordingSolver(solver, num_attr::Integer)
    dim_lyap = length(solver.exponents)
    exponents_history = similar(solver.exponents, (dim_lyap, num_attr))
    LERecordingSolver(solver, exponents_history, 0)
end

@inline function lyapunov_exponents(solver::LERecordingSolver)
    lyapunov_exponents(solver.solver)
end

function step!(solver::LERecordingSolver)
    step!(solver.solver)
    solver.i_hist += 1
    solver.exponents_history[:, solver.i_hist] .= lyapunov_exponents(solver)
end

function solve!(solver::LERecordingSolver; kwargs...)
    _, num_attr = size(solver.exponents_history)
    solver.i_hist = 0
    solve!(solver, num_attr; kwargs...)
end

Base.show(io::IO, solver::LERecordingSolver) = show(io, solver.solver)
