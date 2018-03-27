using DiffEqBase: ODEProblem, RODEProblem


function get_tangent_prob(prob::LEProblem{RODEProblem},
                          u0 = phase_tangent_state(prob);
                          kwargs...)
    phase_prob = prob.phase_prob

    dim_phase = length(phase_prob.u0)
    rand_prototype = if phase_prob.rand_prototype === nothing
        similar(phase_prob.u0, dim_phase)
    else
        phase_prob.rand_prototype
    end

    return remake(phase_prob;
                  f = get_tangent_dynamics(prob),
                  u0 = u0,
                  rand_prototype = rand_prototype,
                  kwargs...)
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
