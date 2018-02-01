function allocate_array_of_arrays(num, inner_size, inner_type,
                                  constructor = inner_type)
    arrays = Array{inner_type, 1}(num)
    for i in 1:num
        arrays[i] = constructor(inner_size)
    end
    return arrays
end

const UTM = UpperTriangular{Float64, Matrix{Float64}}
make_UTM(args...) = UpperTriangular(Matrix{Float64}(args...))
# TODO: UpperTriangular wastes almost half of memory; fix it

mutable struct CLVSolver{Intr, V, M} <: AbstractLESolver{Intr}
    solver::LESolver{Intr, V, M}
    G_history::Vector{M}
    R_history::Vector{UTM}
    C_history::Vector{UTM}
    num_back::Int

    function CLVSolver(solver::LESolver{Intr, V, M},
                       num_attr,
                       ) where {Intr, V, M}
        dim_phase, dim_lyap = size(solver.tangent_state)
        @assert dim_phase == dim_lyap
        dims = (dim_phase, dim_phase)
        G_history = allocate_array_of_arrays(num_attr, dims, M)
        R_history = allocate_array_of_arrays(num_attr, dims, UTM, make_UTM)
        C_history = allocate_array_of_arrays(num_attr, dims, UTM, make_UTM)
        new{Intr, V, M}(
            solver,
            G_history,
            R_history,
            C_history,
            )
    end
end

function solve!(solver::CLVSolver; kwargs...)
    forward!(solver; kwargs...)
    backward!(solver; kwargs...)
end

function forward!(solver::CLVSolver; kwargs...)
    forward!(solver, length(solver.G_history); kwargs...)
end

function step!(solver::CLVSolver)
    step!(solver.solver)
    n = solver.solver.num_orth
    # TODO: these assignments check if lower half triangle is zero; skip that
    solver.G_history[n] .= solver.solver.tangent_state  # aka Q
    solver.R_history[n] .= solver.solver.R
    return solver
end

function backward!(solver::CLVSolver; progress=-1)
    solver.num_back = back_steps = length(solver.G_history)
    eye!(solver.C_history[back_steps])
    @showprogress_if(
        (progress >= 0), progress, "Computing CLVs...",
        for _ in 1:back_steps - 1
            step_back!(solver)
        end)
end

function step_back!(solver::CLVSolver)
    n = (solver.num_back -= 1)
    R = solver.R_history[n]
    C1 = solver.C_history[n + 1]
    C0 = solver.C_history[n]

    # C₀ = R⁻¹ C₁ D  (Eq. 32, Ginelli et al., 2013)
    A_ldiv_B!(R, (C0 .= C1))
    for i in 1:size(C0, 1)
        C0[:, i] /= norm(@view C0[:, i])
    end
end

function clv(solver, n)
    G = solver.G_history[n]
    C = solver.C_history[n]
    return G * C
end
