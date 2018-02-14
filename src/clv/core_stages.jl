# Some "mix-in" methods for generating default `is_finished`:
stage_index(stage::AbstractComputationStage) = stage.i
is_finished(stage::AbstractComputationStage) =
    stage_index(stage) >= stage_length(stage)

record!(::AbstractComputationStage, ::Any) = nothing


mutable struct ForwardRelaxer{T <: LESolver} <: AbstractComputationStage
    le_solver::T
    num_forward_tran::Int
end

ForwardRelaxer(prob::CLVProblem, ::CLVProblem, ::CLVSolution) =
    ForwardRelaxer(prob)
ForwardRelaxer(prob::CLVProblem) = ForwardRelaxer(get_le_solver(prob),
                                                  prob.num_forward_tran)

stage_index(frx::ForwardRelaxer) = frx.le_solver.num_orth
stage_length(frx::ForwardRelaxer) = frx.num_forward_tran
step!(frx::ForwardRelaxer) = step!(frx.le_solver)

mutable struct ForwardDynamics{S <: CLVSolution,
                               T <: LESolver,
                               } <: AbstractComputationStage
    le_solver::T
    i::Int
    sol::S
    R_history::Vector{UTM}

    ForwardDynamics(le_solver::T, sol::S) where {S, T} =
        new{S, T}(le_solver, 0, sol, sol.R_history)
end

ForwardDynamics(frx::ForwardRelaxer, ::CLVProblem, sol::CLVSolution) =
    ForwardDynamics(frx.le_solver, sol)

stage_length(fitr::ForwardDynamics) = length(fitr.R_history)

@inline function step!(fitr::ForwardDynamics)
    step!(fitr.le_solver)
    i = (fitr.i += 1)
    # TODO: this assignment check if lower half triangle is zero; skip that
    fitr.R_history[i] .= CLV.R_prev(fitr)
    record!(fitr, Val{:G})
    record!(fitr, Val{:x})
    return fitr
end

record!(fitr::ForwardDynamics{<:CLVSolG}, ::Type{Val{:G}}) =
    fitr.sol.G_history[fitr.i] .= CLV.G(fitr)

record!(fitr::ForwardDynamics{<:CLVSolX}, ::Type{Val{:x}}) =
    fitr.sol.x_history[fitr.i] .= phase_state(fitr)

const ForwardPass = Union{ForwardRelaxer, ForwardDynamics}

current_result(fitr::ForwardPass) = CLV.G(fitr)
phase_state(fitr::ForwardPass) = phase_state(fitr.le_solver)


mutable struct BackwardRelaxer <: AbstractComputationStage
    num_backward_tran::Int
    R_history::Vector{UTM}
    C::UTM
    i::Int
end

BackwardRelaxer(fitr::ForwardDynamics, prob::CLVProblem, ::CLVSolution) =
    BackwardRelaxer(fitr, prob)
BackwardRelaxer(fitr::ForwardDynamics, prob::CLVProblem) =
    BackwardRelaxer(prob.num_backward_tran,
                    fitr.R_history,
                    UTM(eye(fitr.R_history[end])),  # TODO: randomize
                    -1)

stage_length(brx::BackwardRelaxer) = brx.num_backward_tran


mutable struct BackwardDynamics{with_D,
                                S <: CLVSolution,
                                } <: AbstractComputationStage
    sol::S
    R_history::Vector{UTM}
    C::UTM
    i::Int
    D_diag::Vector{Float64}

    function BackwardDynamics{with_D}(sol::S, R_history, C, i = -1,
                                      ) where {with_D, S}
        if sol isa CLVSolD && ! with_D
            error("with_D (= $with_D) has to be true when :D is",
                  " in record flag.")
        end
        bitr = new{with_D, S}(sol, R_history, C, i)
        if sol isa CLVSolC
            bitr.sol.C_history[end] .= CLV.C(bitr)
        end
        if sol isa CLVSolD
            bitr.sol.D_history[end] .= NaN
        end
        if with_D
            bitr.D_diag = similar(C, size(C, 1))
        end
        return bitr
    end
end

const BackwardDynamicsRecC = BackwardDynamics{with_D, <:CLVSolC} where {with_D}
const BackwardDynamicsRecD = BackwardDynamics{true, <:CLVSolD}
const BackwardDynamicsWithD = BackwardDynamics{true}

BackwardDynamics(brx::BackwardRelaxer, args...) =
    BackwardDynamics{false}(brx, args...)
BackwardDynamics{with_D}(brx::BackwardRelaxer,
                         ::CLVProblem,
                         sol,
                         ) where {with_D} =
    BackwardDynamics{with_D}(
        sol,
        brx.R_history[1:end - brx.i - 1],
        brx.C)

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

    record!(bitr, Val{:C})
    record!(bitr, Val{:D})
end

@inline Dᵢᵢ⁻¹(C::AbstractArray, i) = norm(@view C[:, i])
@inline Dᵢᵢ⁻¹(bitr::BackwardPass, i) = Dᵢᵢ⁻¹(bitr.C, i)
@inline function Dᵢᵢ⁻¹(bitr::BackwardDynamicsWithD, i)
    di = Dᵢᵢ⁻¹(bitr.C, i)
    bitr.D_diag[i] = 1 / di
    return di
end

function record!(bitr::BackwardDynamicsRecC, ::Type{Val{:C}})
    n = length(bitr.sol.C_history) - bitr.i - 1
    if n > 0
        bitr.sol.C_history[n] .= CLV.C(bitr)
    end
    # TODO: maybe save n=0 too
    # TODO: this assignment check if lower half triangle is zero; skip that
end

function record!(bitr::BackwardDynamicsRecD, ::Type{Val{:D}})
    n = length(bitr.sol.D_history) - bitr.i - 1
    if n > 0
        bitr.sol.D_history[n] .= bitr.D_diag
    end
    # TODO: maybe save n=0 too
end

current_result(bitr::BackwardPass) = bitr.C
