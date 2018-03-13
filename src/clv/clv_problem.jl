using DiffEqBase: ODEProblem, DiscreteProblem

"""
    CLVProblem(phase_prob, t_clv; <keyword arguments>)
    CLVProblem(le_prob::LEProblem; <keyword arguments>)

Covariant Lyapunov vector (CLV) problem.  This is a struc that holds
the dynamical system definition (`phase_prob` and `tangent_dynamics!`)
and the configuration parameters for the algorithm (`t_clv`, etc.).

The CLVs are calculated using the 'dynamical' algorithm proposed by
Ginelli et al. (2007, 2013).

### Arguments
- `t_clv::Real`: Time spam for which CLV are sampled.
  It is `0.8 * le_prob.t_attr` when constructed from [`LEProblem`](@ref).
- `t_forward_tran::Real`, `t_backward_tran::Real`:
  Forward and backward transient dynamics.  They are
  `0.1 * le_prob.t_attr` when constructed from [`LEProblem`](@ref).
- See [`LEProblem`](@ref) for `phase_prob`, `t_renorm`,
  `tangent_dynamics!` and `Q0`.

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
struct CLVProblem{DEP,
                  T <: Real,
                  } <: AbstractSource
    phase_prob::DEP
    # t_phase_tran
    t_forward_tran::T
    t_backward_tran::T
    t_clv::T
    t_renorm::T
    Q0
    tangent_dynamics!

    function CLVProblem{DEP}(
            phase_prob::DEP, t_clv;
            t_forward_tran = 0,
            t_backward_tran = 0,
            t_renorm = 1,
            dim_lyap::Int = dimension(phase_prob),
            Q0 = default_Q0(phase_prob, dimension(phase_prob), dim_lyap),
            tangent_dynamics! = nothing,
            ) where {DEP}

        if t_clv / t_renorm < 1
            error("t_clv / t_renorm = $(t_clv / t_renorm)",
                  " (given: t_clv=$t_clv, t_renorm=$t_renorm)",
                  " must be larger than 1.")
        end

        dim_phase = dimension(phase_prob)
        Q0_size_1 = size(Q0, 1)
        if Q0_size_1 != dim_phase
            error("size(Q0, 1) = $Q0_size_1 != dim_phase = $dim_phase")
        end

        if ! is_semi_unitary(Q0)
            error("Columns in Q0 are not orthonormal.")
        end


        TT = time_type((DEP <: DiscreteProblem);
                       t_forward_tran = t_forward_tran,
                       t_backward_tran = t_backward_tran,
                       t_clv = t_clv,
                       t_renorm = t_renorm)

        new{DEP, TT}(
            phase_prob,
            t_forward_tran, t_backward_tran, t_clv, t_renorm,
            Q0, tangent_dynamics!)
    end
end

CLVProblem(phase_prob::DEP, args...; kwargs...) where {DEP <: ODEProblem} =
    CLVProblem{ODEProblem}(phase_prob, args...; kwargs...)
CLVProblem(phase_prob::DEP, args...; kwargs...) where {DEP <:DiscreteProblem} =
    CLVProblem{DiscreteProblem}(phase_prob, args...; kwargs...)

CLVProblem(prob::LEProblem{DEP};
           t_clv = ceil_if(DEP <: DiscreteProblem, prob.t_attr * 0.8),
           t_renorm = prob.t_renorm,
           kwargs...) where {DEP} =
    CLVProblem(
        prob.phase_prob, t_clv;
        # t_phase_tran = prob.t_tran,
        t_forward_tran = ceil_if(DEP <: DiscreteProblem, prob.t_attr * 0.1),
        t_backward_tran = ceil_if(DEP <: DiscreteProblem, prob.t_attr * 0.1),
        t_renorm = t_renorm,
        Q0 = prob.Q0,
        tangent_dynamics! = prob.tangent_dynamics!,
        kwargs...)

# DRY: it's almost the same as for prob::LEProblem
function phase_tangent_state(prob::CLVProblem, x0 = prob.phase_prob.u0)
    dim_phase = dimension(prob.phase_prob)
    Q0 = prob.Q0
    dim_phase, dim_lyap = size(Q0)
    u0 = similar(x0, (dim_phase, dim_lyap + 1))
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
                       t_renorm = prob.t_renorm,
                       kwargs...)
    tangent_prob = get_tangent_prob(prob, u0)
    t_attr = prob.t_forward_tran + prob.t_clv + prob.t_backward_tran
    return TangentRenormalizer(get_integrator(tangent_prob),
                               t_attr, t_renorm; kwargs...)
end


mutable struct CLVSolution{RecG, RecC, RecD, RecX}
    num_clv::Int
    num_backward_tran::Int
    # TODO: make those types more general
    R_history::Vector{UTM}
    G_history::Vector{Matrix{Float64}}
    C_history::Vector{UTM}
    D_history::Vector{Vector{Float64}}
    x_history::Vector{Vector{Float64}}

    # TODO: Do not store backward transient part in CLVSolution.
    # Separate it out to some kind of cache object.

    function CLVSolution{RecG, RecC, RecD, RecX}(
            dim_phase::Int, dim_lyap::Int,
            num_clv::Int, num_backward_tran::Int,
            ) where {RecG, RecC, RecD, RecX}

        sol = new{RecG, RecC, RecD, RecX}(num_clv, num_backward_tran)

        sol.R_history = allocate_array_of_arrays(
            num_clv + num_backward_tran,
            (dim_lyap, dim_lyap),
            UTM, make_UTM)
        # Maybe it makes sense to have an option to progressively
        # delete R during the backward iteration to save some memory.

        if RecG
            GT = Matrix{Float64}
            sol.G_history = allocate_array_of_arrays(
                num_clv + num_backward_tran,  # TOOD: remove num_backward_tran
                (dim_phase, dim_lyap),
                GT)
        end

        if RecC
            sol.C_history = allocate_array_of_arrays(
                num_clv,
                (dim_lyap, dim_lyap),
                UTM, make_UTM)
        end

        if RecD
            DT = Vector{Float64}
            sol.D_history = allocate_array_of_arrays(
                num_clv,
                (dim_lyap,),
                DT)
        end

        if RecX
            XT = Vector{Float64}
            sol.x_history = allocate_array_of_arrays(
                num_clv + num_backward_tran,  # TOOD: remove num_backward_tran
                (dim_phase,),
                XT)
        end

        # TODO: Do not save G and x for num_backward_tran.
        # Separate ForwardDynamics into two parts.

        return sol
    end
end

function CLVSolution(prob::CLVProblem, record::Vector{Symbol})
    unsupported = setdiff(record, [:G, :C, :D, :x])
    if ! isempty(unsupported)
        error("Unsupported record key(s): $unsupported")
    end

    dim_phase, dim_lyap = size(prob.Q0)
    return CLVSolution{
        (:G in record),
        (:C in record),
        (:D in record),
        (:x in record),
    }(dim_phase,
      dim_lyap,
      ceil(Int, prob.t_clv / prob.t_renorm),
      ceil(Int, prob.t_backward_tran / prob.t_renorm))
end

# const CLVSolR = CLVSolution
const CLVSolG = CLVSolution{true}
const CLVSolC = CLVSolution{RecG, true} where {RecG}
const CLVSolD = CLVSolution{RecG, RecC, true} where {RecG, RecC}
const CLVSolX = CLVSolution{RecG, RecC, RecD, true} where {RecG, RecC, RecD}
