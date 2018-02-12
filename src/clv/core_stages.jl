# Some "mix-in" methods for generating default `is_finished`:
function stage_length end
stage_index(stage::AbstractComputationStage) = stage.i
is_finished(stage::AbstractComputationStage) =
    stage_index(stage) >= stage_length(stage)

record!(::AbstractComputationStage, ::Any) = nothing
# TODO: Separate recording code into "aggregator" types and move
# flag-based dispatching there.  This way, stages have to depend on
# single type parameter and code noise here would be reduced.


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

mutable struct ForwardDynamics{record_G,
                               record_x,
                               T <: LESolver,
                               } <: AbstractComputationStage
    le_solver::T
    i::Int
    sol::CLVSolution
    R_history::Vector{UTM}

    ForwardDynamics{record_G,
                    record_x,
                    }(le_solver::T, i, sol) where {record_G, record_x, T} =
        new{record_G, record_x, T}(le_solver, i, sol)
end

ForwardDynamics_from_record(record) = ForwardDynamics{(:G in record),
                                                      (:x in record)}
ForwardDynamics(a, b, c) = ForwardDynamics{false, false}(a, b, c)
ForwardDynamicsWithGHistory(a, b, c) = ForwardDynamics{true}(a, b, c)

ForwardDynamics{record_G, record_x}(frx::ForwardRelaxer,
                                    prob::CLVProblem, sol,
                                    ) where {record_G, record_x} =
    ForwardDynamics{record_G, record_x}(frx.le_solver, prob::CLVProblem, sol)
function ForwardDynamics{record_G,
                         record_x}(le_solver::LESolver, prob::CLVProblem, sol,
                                   ) where {record_G, record_x}
    fitr = ForwardDynamics{record_G, record_x}(le_solver, 0, sol)
    allocate_solution!(fitr, prob)
    return fitr
end

function allocate_solution!(fitr::ForwardDynamics{record_G,
                                                  record_x},
                            prob) where {record_G, record_x}

    dim_phase, dim_lyap = size(fitr.le_solver.tangent_state)
    fitr.sol.R_history = allocate_array_of_arrays(
        prob.num_clv + prob.num_backward_tran,
        (dim_lyap, dim_lyap),
        UTM, make_UTM)
    fitr.R_history = fitr.sol.R_history

    if record_G
        # When CLVSolution is parameterized, GT has to be extracted
        # from there:
        GT = Matrix{Float64}
        fitr.sol.G_history = allocate_array_of_arrays(
            prob.num_clv + prob.num_backward_tran,
            size(fitr.le_solver.tangent_state),
            GT)
    end

    if record_x
        XT = Vector{Float64}
        fitr.sol.x_history = allocate_array_of_arrays(
            prob.num_clv + prob.num_backward_tran,
            size(fitr.le_solver.phase_state),
            XT)
    end
    # TODO: Do not save G and x for prob.num_backward_tran; Separate
    # ForwardDynamics into two parts.
end

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

record!(fitr::ForwardDynamics{true}, ::Type{Val{:G}}) =
    fitr.sol.G_history[fitr.i] .= CLV.G(fitr)

record!(fitr::ForwardDynamics{record_G, true},
        ::Type{Val{:x}}) where {record_G} =
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
                                record_C,
                                } <: AbstractComputationStage
    sol::CLVSolution
    R_history::Vector{UTM}
    C::UTM
    i::Int
    D_diag::Vector{Float64}

    function BackwardDynamics{with_D, record_C,
                              }(sol, R_history, C, i = -1,
                                ) where {with_D, record_C}
        bitr = new{with_D, record_C}(sol, R_history, C, i)
        if with_D
            bitr.D_diag = similar(C, size(C, 1))
        end
        return bitr
    end
end

BackwardDynamics_from_record(record) = BackwardDynamics{false, (:C in record)}

const BackwardDynamicsWithD = BackwardDynamics{true, false}
const BackwardDynamicsWithCHistory = BackwardDynamics{false, true}

BackwardDynamics(brx::BackwardRelaxer, args...) =
    BackwardDynamics{false, false}(brx, args...)
BackwardDynamics{with_D, record_C}(brx::BackwardRelaxer,
                                   ::CLVProblem,
                                   sol,
                                   ) where {with_D, record_C} =
    BackwardDynamics{with_D, record_C}(brx, sol)

function BackwardDynamics{with_D, record_C}(brx::BackwardRelaxer,
                                            sol,
                                            ) where {with_D, record_C}
    bitr = BackwardDynamics{with_D, record_C}(
        sol,
        brx.R_history[1:end - brx.i - 1],
        brx.C)
    allocate_solution!(bitr)
    return bitr
end

function allocate_solution!(bitr::BackwardDynamics{with_D, record_C},
                            ) where {with_D, record_C}
    if record_C
        bitr.sol.C_history = allocate_array_of_arrays(length(bitr),
                                                      size(bitr.C),
                                                      UTM, make_UTM)
        bitr.sol.C_history[end] = CLV.C(bitr)
    end
end

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
end

function record!(bitr::BackwardDynamics{with_D, true},
                 ::Type{Val{:C}}) where {with_D}
    n = length(bitr.sol.C_history) - bitr.i - 1
    if n > 0
        bitr.sol.C_history[n] .= CLV.C(bitr)
    end
    # TODO: maybe save n=0 too
    # TODO: this assignment check if lower half triangle is zero; skip that
end

@inline Dᵢᵢ⁻¹(C::AbstractArray, i) = norm(@view C[:, i])
@inline Dᵢᵢ⁻¹(bitr::BackwardPass, i) = Dᵢᵢ⁻¹(bitr.C, i)
@inline function Dᵢᵢ⁻¹(bitr::BackwardDynamics{true}, i)
    di = Dᵢᵢ⁻¹(bitr.C, i)
    bitr.D_diag[i] = 1 / di
    return di
end

current_result(bitr::BackwardPass) = bitr.C
