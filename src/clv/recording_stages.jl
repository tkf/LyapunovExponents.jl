abstract type AbstractRecordingStage <: AbstractComputationStage end

is_finished(stage::AbstractRecordingStage) = is_finished(stage.stage)
# stage_index(stage::AbstractRecordingStage) = stage_index(stage.stage)
# stage_length(stage::AbstractRecordingStage) = stage_length(stage.stage)

struct ForwardDynamicsWithGHistory <: AbstractRecordingStage
    stage::ForwardDynamics
    G_history::Vector{Matrix{Float64}}  # TODO: make it more general
end

ForwardDynamicsWithGHistory(frx::ForwardRelaxer, prob::CLVProblem) =
    ForwardDynamicsWithGHistory(
        ForwardDynamics(frx, prob),
        allocate_forward_history(frx.le_solver, prob, Matrix{Float64}),
    )

@inline function step!(fitr::ForwardDynamicsWithGHistory)
    step!(fitr.stage)
    # TODO: this assignment check if lower half triangle is zero; skip that
    fitr.G_history[stage_index(fitr.stage)] .= CLV.G(fitr)
    # This is Q of the QR decomposition.
end

record!(sol::CLVSolution, fitr::ForwardDynamicsWithGHistory) =
    sol.G_history = fitr.G_history

BackwardRelaxer(fitr::ForwardDynamicsWithGHistory, prob::CLVProblem) =
    BackwardRelaxer(fitr.stage, prob)
# TODO: Find a way (trait?) to avoid writing this for all combinations
# of different stages.


# struct ForwardDynamicsWithCheckPoints <: AbstractComputationStage
#     iter::ForwardDynamics
#     ???
# end


struct BackwardDynamicsWithCHistory <: AbstractRecordingStage
    stage::BackwardDynamics
    C_history::Vector{UTM}
end

BackwardDynamicsWithCHistory(brx::BackwardRelaxer, _) =
    BackwardDynamicsWithCHistory(brx)
BackwardDynamicsWithCHistory(brx::BackwardRelaxer) =
    BackwardDynamicsWithCHistory(
        BackwardDynamics(brx),
        allocate_array_of_arrays(length(brx.R_history) - brx.i,
                                 size(brx.C),
                                 UTM, make_UTM),
    )

@inline function step!(bitr::BackwardDynamicsWithCHistory)
    step!(bitr.stage)
    # TODO: this assignment check if lower half triangle is zero; skip that
    bitr.C_history[stage_index(bitr.stage)] .= CLV.C(bitr)
end

record!(sol::CLVSolution, bitr::BackwardDynamicsWithCHistory) =
    sol.C_history = bitr.C_history
