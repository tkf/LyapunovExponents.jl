using DiffEqBase: ODEProblem, RODEProblem

function make_tangent_prob(phase_prob::RODEProblem, tangent_dynamics, u0)
    dim_phase = length(phase_prob.u0)
    rand_prototype = if phase_prob.rand_prototype === nothing
        similar(phase_prob.u0, dim_phase)
    else
        phase_prob.rand_prototype
    end

    return remake(
        phase_prob;
        f = tangent_dynamics,
        u0 = u0,
        rand_prototype = rand_prototype)
end


function with_output_noise(ode::ODEProblem,
                           g::Union{AbstractVector, Number};
                           kwargs...)
    function f(du, u, p, t, W)
        ode.f(du, u, p, t)
        du .+= g .* W
    end

    return RODEProblem(f, ode.u0, ode.tspan, ode.p; kwargs...)
end


function with_input_noise(ode::ODEProblem,
                          g::Union{AbstractVector, Number};
                          kwargs...)
    function f(du, u, p, t, W)
        u .+= g .* W
        ode.f(du, u, p, t)
    end

    return RODEProblem(f, ode.u0, ode.tspan, ode.p; kwargs...)
end
