get_dim_lyap(integrator::Union{DEIntegrator, DiscreteIterator}) =
    size(init_tangent_state(integrator))[2]

function get_dim_lyap(prob::Union{LEProblem, CLVProblem})
    Q0 = @view prob.tangent_prob.u0[:, 2:end]
    return size(Q0, 2)
end
