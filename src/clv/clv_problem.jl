struct CLVProblem <: AbstractStage
end

finish!(::CLVProblem) = nothing
is_finished(::CLVProblem) = true

CLVProblem(prob::LEProblem) =
    CLVProblem()
