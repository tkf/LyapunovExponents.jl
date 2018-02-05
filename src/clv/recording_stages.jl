abstract type AbstractRecordingStage <: AbstractComputationStage end

is_finished(stage::AbstractRecordingStage) = is_finished(stage.stage)
# stage_index(stage::AbstractRecordingStage) = stage_index(stage.stage)
# stage_length(stage::AbstractRecordingStage) = stage_length(stage.stage)

struct ForwardDynamicsWithGHistory <: AbstractRecordingStage
    stage::ForwardDynamics
    G_history::Vector{Matrix{Float64}}  # TODO: make it more general
end

@inline function step!(fitr::ForwardDynamicsWithGHistory)
    step!(fitr.stage)
    # TODO: this assignment check if lower half triangle is zero; skip that
    G_history[stage_index(fitr.stage)] .= CLV.G(fitr)
    # This is Q of the QR decomposition.
end

record!(sol::CLVSolution, fitr::ForwardDynamicsWithGHistory) =
    sol.G_history = fitr.G_history


# struct ForwardDynamicsWithCheckPoints <: AbstractComputationStage
#     iter::ForwardDynamics
#     ???
# end


struct BackwardDynamicsWithCHistory <: AbstractRecordingStage
    bitr::BackwardDynamics
    C_history::Vector{UTM}
end

@inline function step!(bitr::BackwardDynamicsWithCHistory)
    step!(bitr.stage)
    # TODO: this assignment check if lower half triangle is zero; skip that
    bitr.C_history[stage_index(bitr.stage)] .= CLV.C(bitr)
end
