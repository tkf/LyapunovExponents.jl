# Some "mix-in" methods for generating default `is_finished`:
function stage_length end
stage_index(stage::AbstractComputationStage) = stage.i
is_finished(stage::AbstractComputationStage) =
    stage_index(stage) >= stage_length(stage)


mutable struct ForwardRelaxer{T <: LESolver} <: AbstractComputationStage
    le_solver::T
    num_forward_tran::Int
end

ForwardRelaxer(prob::CLVProblem, _) = ForwardRelaxer(prob)
ForwardRelaxer(prob::CLVProblem) = ForwardRelaxer(get_le_solver(prob),
                                                  prob.num_forward_tran)

stage_index(frx::ForwardRelaxer) = frx.le_solver.num_orth
stage_length(frx::ForwardRelaxer) = frx.num_forward_tran
step!(frx::ForwardRelaxer) = step!(frx.le_solver)

mutable struct ForwardDynamics{T <: LESolver} <: AbstractComputationStage
    le_solver::T
    R_history::Vector{UTM}
    i::Int
end

ForwardDynamics(frx::ForwardRelaxer, prob::CLVProblem) =
    ForwardDynamics(frx.le_solver, prob::CLVProblem)
ForwardDynamics(le_solver::LESolver, prob::CLVProblem) =
    ForwardDynamics(le_solver,
                    allocate_forward_history(le_solver, prob, UTM, make_UTM),
                    0)

function allocate_forward_history(le_solver::LESolver, prob::CLVProblem,
                                  inner_type, args...)
    num = prob.num_clv + prob.num_backward_tran
    dim_phase, dim_lyap = size(le_solver.tangent_state)
    @assert dim_phase == dim_lyap
    dims = (dim_phase, dim_phase)
    return allocate_array_of_arrays(num, dims, inner_type, args...)
end

stage_length(fitr::ForwardDynamics) = length(fitr.R_history)

@inline function step!(fitr::ForwardDynamics)
    step!(fitr.le_solver)
    i = (fitr.i += 1)
    # TODO: this assignment check if lower half triangle is zero; skip that
    fitr.R_history[i] .= CLV.R_prev(fitr)
    return fitr
end

const ForwardPass = Union{ForwardRelaxer, ForwardDynamics}

current_result(fitr::ForwardPass) = CLV.R_prev(fitr)
phase_state(fitr::ForwardPass) = phase_state(fitr.le_solver)


mutable struct BackwardRelaxer <: AbstractComputationStage
    num_backward_tran::Int
    R_history::Vector{UTM}
    C::UTM
    i::Int
end

BackwardRelaxer(fitr::ForwardDynamics, prob::CLVProblem) =
    BackwardRelaxer(prob.num_backward_tran,
                    fitr.R_history,
                    UTM(eye(fitr.R_history[end])),  # TODO: randomize
                    -1)

stage_length(brx::BackwardRelaxer) = brx.num_backward_tran


mutable struct BackwardDynamics{with_D} <: AbstractComputationStage
    R_history::Vector{UTM}
    C::UTM
    i::Int
    D_diag::Vector{Float64}

    function BackwardDynamics{with_D}(R_history, C, i = -1) where {with_D}
        bitr = new{with_D}(R_history, C, i)
        if with_D
            bitr.D_diag = similar(C, size(C, 1))
        end
        return bitr
    end
end

const BackwardDynamicsWithD = BackwardDynamics{true}

BackwardDynamics(brx::BackwardRelaxer, args...) =
    BackwardDynamics{false}(brx, args...)
BackwardDynamics{with_D}(brx::BackwardRelaxer, _) where {with_D} =
    BackwardDynamics{with_D}(brx)
BackwardDynamics{with_D}(brx::BackwardRelaxer) where {with_D} =
    BackwardDynamics{with_D}(brx.R_history[1:end - brx.i - 1], brx.C)

stage_length(bitr::BackwardDynamics) = length(bitr.R_history)

const BackwardPass = Union{BackwardRelaxer, BackwardDynamics}

stage_index(stage::BackwardPass) = stage.i + 1  # TODO: don't

@inline function step!(bitr::BackwardPass)
    i = (bitr.i += 1)
    R = bitr.R_history[end-i]
    C = bitr.C

    # C₀ = R⁻¹ C₁ D  (Eq. 32, Ginelli et al., 2013)
    A_ldiv_B!(R, C)
    # now:  C = R⁻¹ C₁
    for i in 1:size(C, 1)
        # C[:, i] /= norm(@view C[:, i])
        C[:, i] /= Dᵢᵢ⁻¹(bitr, i)
    end
    # now:  C = C₀ = R⁻¹ C₁ D
end

@inline Dᵢᵢ⁻¹(C::AbstractArray, i) = norm(@view C[:, i])
@inline Dᵢᵢ⁻¹(bitr::BackwardPass, i) = Dᵢᵢ⁻¹(bitr.C, i)
@inline function Dᵢᵢ⁻¹(bitr::BackwardDynamicsWithD, i)
    di = Dᵢᵢ⁻¹(bitr.C, i)
    bitr.D_diag[i] = 1 / di
    return di
end

current_result(bitr::BackwardPass) = bitr.C
