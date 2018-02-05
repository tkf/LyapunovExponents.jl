# Some "mix-in" methods for generating default `is_finished`:
function stage_length end
stage_index(stage::AbstractComputationStage) = stage.i
is_finished(stage::AbstractComputationStage) =
    stage_index(stage) + 1 >= stage_length(stage)


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
    fitr.R_history[end-i] .= CLV.R(fitr)
    return fitr
end

const ForwardPass = Union{ForwardRelaxer, ForwardDynamics}

current_result(fitr::ForwardPass) = CLV.R(fitr)
phase_state(fitr::ForwardPass) = phase_state(fitr.le_solver)


mutable struct BackwardRelaxer <: AbstractComputationStage
    num_backward_tran
    R_history
    C
    i
end

BackwardRelaxer(fitr::ForwardDynamics, prob::CLVProblem) =
    BackwardRelaxer(prob.num_backward_tran,
                    fitr.R_history,
                    eye(fitr.R_history[end]),  # TODO: randomize
                    0)

stage_length(brx::BackwardRelaxer) = brx.num_backward_tran


mutable struct BackwardDynamics <: AbstractComputationStage
    R_history
    C
    i
end

BackwardDynamics(brx::BackwardRelaxer, _) = BackwardDynamics(brx)
BackwardDynamics(brx::BackwardRelaxer) =
    BackwardDynamics(brx.R_history[1:end - brx.i], brx.C, 0)

stage_length(bitr::BackwardDynamics) = length(bitr.R_history)

const BackwardPass = Union{BackwardRelaxer, BackwardDynamics}

@inline function step!(bitr::BackwardPass)
    i = (bitr.i += 1)
    R = bitr.R_history[end-i]
    C = bitr.C

    # C₀ = R⁻¹ C₁ D  (Eq. 32, Ginelli et al., 2013)
    A_ldiv_B!(R, C)
    for i in 1:size(C, 1)
        C[:, i] /= norm(@view C[:, i])
    end
end

current_result(bitr::BackwardPass) = bitr.C
