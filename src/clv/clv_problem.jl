using DiffEqBase: DEProblem, ODEProblem, DiscreteProblem, parameterless_type

"""
    CLVProblem(phase_prob, t_clv; <keyword arguments>)
    CLVProblem(le_prob::LEProblem; <keyword arguments>)

Covariant Lyapunov vector (CLV) problem.  This is a struct that holds
the dynamical system definition (`phase_prob` and `tangent_dynamics`)
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
  `tangent_dynamics` and `Q0`.

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
struct CLVProblem{PPr, TPr,
                  T <: Real,
                  } <: AbstractSource
    phase_prob::PPr
    tangent_prob::TPr
    # t_phase_tran
    t_forward_tran::T
    t_backward_tran::T
    t_clv::T
    t_renorm::T

    function CLVProblem{PPr, TPr}(
            phase_prob::PPr,
            tangent_prob::TPr,
            t_clv;
            t_forward_tran = 0,
            t_backward_tran = 0,
            t_renorm = 1,
            ) where {PPr, TPr}

        if t_clv / t_renorm < 1
            error("t_clv / t_renorm = $(t_clv / t_renorm)",
                  " (given: t_clv=$t_clv, t_renorm=$t_renorm)",
                  " must be larger than 1.")
        end

        TT = time_type((PPr <: DiscreteProblem);
                       t_forward_tran = t_forward_tran,
                       t_backward_tran = t_backward_tran,
                       t_clv = t_clv,
                       t_renorm = t_renorm)

        new{PPr, TPr, TT}(
            phase_prob, tangent_prob,
            t_forward_tran, t_backward_tran, t_clv, t_renorm)
    end
end

function CLVProblem(phase_prob::PPr;
                    t_clv::Real = error("Keyword argument",
                                        " `t_clv` is required."),
                    tangent_dynamics = nothing,
                    tangent_prob::Union{Void, DEProblem} = nothing,
                    dim_lyap = dimension(phase_prob),
                    Q0 = default_Q0(phase_prob, dimension(phase_prob), dim_lyap),
                    kwargs...) where {PPr <: DEProblem}
    tangent_prob = validate_tangent_prob(
        phase_prob, tangent_prob, tangent_dynamics, Q0)
    CLVProblem{parameterless_type(PPr),
              parameterless_type(typeof(tangent_prob))}(
        phase_prob, tangent_prob, t_clv;
        kwargs...)
end

CLVProblem(prob::LEProblem{DEP}; kwargs...) where {DEP} =
    CLVProblem(
        prob.phase_prob;
        t_clv = ceil_if(DEP <: DiscreteProblem, prob.t_attr * 0.8),
        # t_phase_tran = prob.t_tran,
        t_forward_tran = ceil_if(DEP <: DiscreteProblem, prob.t_attr * 0.1),
        t_backward_tran = ceil_if(DEP <: DiscreteProblem, prob.t_attr * 0.1),
        t_renorm = prob.t_renorm,
        tangent_prob = prob.tangent_prob,
        kwargs...)

function get_le_solver(prob;
                       x0 = prob.phase_prob.u0,
                       t_renorm = prob.t_renorm,
                       kwargs...)

    tspan = if prob.tangent_prob isa ODEProblem
        (prob.tangent_prob.tspan[1], Inf)
    else
        prob.phase_prob.tspan
    end
    # See: [[../core.jl::get_tangent_integrator]]

    tangent_prob = get_tangent_prob(prob; x0=x0, tspan=tspan)
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

    Q0 = @view prob.tangent_prob.u0[:, 2:end]
    dim_phase, dim_lyap = size(Q0)
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
