mutable struct LERecordingSolver{Intr} <: AbstractLESolver{Intr}
    solver
    exponents_history
    i_hist
end

function LERecordingSolver(solver::AbstractLESolver{Intr},
                           num_attr::Integer) where {Intr}
    exps = solver.inst_exponents
    dim_lyap = length(exps)
    exponents_history = similar(exps, (dim_lyap, num_attr))
    LERecordingSolver{Intr}(solver, exponents_history, 0)
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
