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


"""
    CLVSolver(solver::LESolver, num_attr::Int)
    CLVSolver(demo::LEDemo)

A covariant Lyapunov vector (CLV) calculator.

This solver implements the 'dynamical' algorithm proposed by Ginelli
et al. (2007, 2013).

### Example

Figure 1 of Ginelli et al. (2007) can be reproduced as follows.  This
code first calculates the CLV of the Henon map and plot the histogram
of angle between two CLVs:

```julia
using LyapunovExponents

demo = LyapunovExponents.henon_map(num_attr=100000)
solver = LyapunovExponents.CLVSolver(demo)
@time solve!(solver)

limit = floor(Int, length(solver.C_history) * 0.9)
angles = [acos(abs(dot(C[:, 1], C[:, 2]))) * 2 / π
          for C in solver.C_history[1:limit]]

using Plots
stephist(angles, bins=1000, normalize=true,
         xlabel="Angle [pi/2 rad]", ylabel="Density", label="")
```

### Reference

* Ginelli, F., Poggi, P., Turchi, A., Chaté, H., Livi, R., & Politi, A. (2007).
  *Characterizing Dynamics with Covariant Lyapunov Vectors.*
  Physical Review Letters, 99(13), 1–4.
  <http://doi.org/10.1103/PhysRevLett.99.130601>
* Ginelli, F., Chaté, H., Livi, R., & Politi, A. (2013).
  *Covariant Lyapunov vectors.*
  Journal of Physics A: Mathematical and Theoretical, 46(25), 254005.
  <http://doi.org/10.1088/1751-8113/46/25/254005>
"""
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
