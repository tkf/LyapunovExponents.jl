mutable struct LERecordingSolver{Intr} <: AbstractLESolver{Intr}
    solver
    ftle_history
    i_hist
end

function LERecordingSolver(solver::AbstractLESolver{Intr}) where {Intr}
    num_attr = solver.num_attr
    exps = ftle(solver)
    dim_lyap = length(exps)
    ftle_history = similar(exps, (dim_lyap, num_attr))
    LERecordingSolver{Intr}(solver, ftle_history, 0)
end

@inline function lyapunov_exponents(solver::LERecordingSolver)
    lyapunov_exponents(solver.solver)
end

function step!(solver::LERecordingSolver)
    step!(solver.solver)
    solver.i_hist += 1
    solver.ftle_history[:, solver.i_hist] .= ftle(solver.solver)
end

function solve!(solver::LERecordingSolver; kwargs...)
    _, num_attr = size(solver.ftle_history)
    solver.i_hist = 0
    forward!(solver, num_attr; kwargs...)
end

function exponents_history(solver::LERecordingSolver)
    dim_lyap, num_attr = size(solver.ftle_history)
    ftle_history = solver.ftle_history

    m = VecMean(dim_lyap)
    s = OnlineStats.Series(m)
    le_hist = similar(ftle_history)
    for i in 1:num_attr
        OnlineStats.fit!(s, @view ftle_history[:, i])
        le_hist[:, i] .= mean(m)
    end
    return le_hist
end

Base.show(io::IO, solver::LERecordingSolver) = show(io, solver.solver)
