using DiffEqBase: ODEProblem, DiscreteProblem

"""
    CLVProblem(phase_prob, num_clv; <keyword arguments>)
    CLVProblem(le_prob::LEProblem; <keyword arguments>)

Covariant Lyapunov vector (CLV) problem.  This is a struc that holds
the dynamical system definition (`phase_prob` and `tangent_dynamics!`)
and the configuration parameters for the algorithm (`num_clv`, etc.).

The CLVs are calculated using the 'dynamical' algorithm proposed by
Ginelli et al. (2007, 2013).

### Arguments
- `num_clv::Int`: Number of points at which CLV are sampled.
  It is `0.8 * le_prob.num_attr` when constructed from [`LEProblem`](@ref).
- `num_forward_tran::Int`, `num_backward_tran::Int`:
  Forward and backward transient dynamics.  They are
  `0.1 * le_prob.num_attr` when constructed from [`LEProblem`](@ref).
- See [`LEProblem`](@ref) for `phase_prob`, `tangent_dynamics!` and `Q0`.

### Examples (in the online documentation)
- [Ginelli et al. (2007), Figure 1a](@ref)
- [Covariant Lyapunov vectors on the Lorenz attractor](@ref)

See: <https://tkf.github.io/LyapunovExponents.jl/latest/gallery/>

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
struct CLVProblem{DEP} <: AbstractStage
    phase_prob::DEP
    # num_phase_tran
    num_forward_tran::Int
    num_backward_tran::Int
    num_clv::Int
    Q0
    tangent_dynamics!

    function CLVProblem{DEP}(
            phase_prob::DEP, num_clv::Int;
            num_forward_tran::Int = 0,
            num_backward_tran::Int = 0,
            Q0 = default_Q0(phase_prob, dimension(phase_prob)),
            tangent_dynamics! = nothing,
            ) where {DEP}

        if num_clv < 1
            error("num_clv (given: $num_clv) must be larger than 1.")
        end

        dim_phase = dimension(phase_prob)
        Q0_size = (dim_phase, dim_phase)
        if size(Q0) != Q0_size
            error("Q0 must be a square matrix of size $Q0_size")
        end

        if ! is_semi_unitary(Q0)
            error("Columns in Q0 are not orthonormal.")
        end

        new{DEP}(phase_prob,
                 num_forward_tran, num_backward_tran, num_clv,
                 Q0, tangent_dynamics!)
    end
end

CLVProblem(phase_prob::DEP, num_clv; kwargs...) where {DEP <: ODEProblem} =
    CLVProblem{ODEProblem}(phase_prob, num_clv; kwargs...)
CLVProblem(phase_prob::DEP, num_clv; kwargs...) where {DEP <:DiscreteProblem} =
    CLVProblem{DiscreteProblem}(phase_prob, num_clv; kwargs...)

finish!(::CLVProblem) = nothing
is_finished(::CLVProblem) = true

CLVProblem(prob::LEProblem;
           num_clv = max(1, floor(Int, prob.num_attr * 0.8)),
           kwargs...) =
    CLVProblem(
        prob.phase_prob, num_clv;
        # num_phase_tran = prob.num_tran,
        num_forward_tran = floor(Int, prob.num_attr * 0.1),
        num_backward_tran = floor(Int, prob.num_attr * 0.1),
        Q0 = prob.Q0,
        tangent_dynamics! = prob.tangent_dynamics!,
        kwargs...)

# DRY: it's almost the same as for prob::LEProblem
function phase_tangent_state(prob::CLVProblem, x0 = prob.phase_prob.u0)
    dim_phase = dimension(prob.phase_prob)
    Q0 = prob.Q0
    u0 = similar(x0, (length(x0), dim_phase + 1))
    u0[:, 1] = x0
    u0[:, 2:end] = Q0
    u0
end

# DRY: it's exactly the same as for prob::LEProblem
function get_tangent_prob(prob::CLVProblem{DEP},
                          u0 = phase_tangent_state(prob)) where {DEP}
    phase_prob = prob.phase_prob
    return DEP(
        get_tangent_dynamics(prob),
        u0,
        phase_prob.tspan,
        phase_prob.p,
    )
end

function get_le_solver(prob, u0 = phase_tangent_state(prob);
                       kwargs...)
    tangent_prob = get_tangent_prob(prob, u0)
    num_attr = prob.num_forward_tran + prob.num_clv + prob.num_backward_tran
    return get_le_solver(tangent_prob, num_attr; kwargs...)
end


mutable struct CLVSolution
    # TODO: make those types more general
    G_history::Vector{Matrix{Float64}}
    R_history::Vector{UTM}
    C_history::Vector{UTM}
    x_history::Vector{Vector{Float64}}

    CLVSolution() = new()
end
