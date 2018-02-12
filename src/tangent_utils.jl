# Like set_u0, but more limited and usable for DiscreteProblem
with_u0(prob::ODEProblem, u0) = ODEProblem(prob.f, u0, prob.tspan, prob.p)
with_u0(prob::DiscreteProblem, u0) =
    DiscreteProblem(prob.f, u0, prob.tspan, prob.p)

function augmented_vector(x0::AbstractVector, Q0::AbstractArray)
    u0 = similar(x0, (size(x0, 1), size(Q0, 2) + 1))
    u0[:, 1] = x0
    u0[:, 2:end] = Q0
    return u0
end

augmented_vector(x0::AbstractVector, Q0::UniformScaling) =
    augmented_vector(x0, Q0 * eye(size(x0, 1)))

doc"""
    tangent_propagate(prob_or_solver::Union{DEProblem,LESolver}
                      [, tangent_state];
                      phase_state,
                      <keyword arguments>)

Propagate `tangent_state` according to tangent dynamics defined by
`prob_or_solver`.  Supplying `I` to `tangent_state` is a quick (but
not efficient) way to obtain the cocycle ``M_{k,n}`` (tangent space
"evolution operator") which is:

* for the maps: ``M_{k,n} = \prod_{i=0}^{i<k} df/dx(n+i)``.
* for the flows: the solution ``M_{k,n} = M(t_{k+n})`` of the linear
  ODE ``dM/dt = (df/dx(t)) M`` with the initial condition
  ``M(t_{n}) = I``.
"""
tangent_propagate(prob_or_solver::Union{DEProblem,LESolver},
                  args...; kwargs...) =
    tangent_propagate(Val{:tangent}, prob_or_solver, args...; kwargs...)

tangent_propagate(::Type{Val{:tangent}}, args...; kwargs...) =
    @view tangent_propagate(Val{:full}, args...; kwargs...)[:, 2:end]

tangent_propagate(::Type{Val{:full}}, args...; kwargs...) =
    current_state(tangent_propagate(Val{:integrator}, args...; kwargs...))

function tangent_propagate(::Type{Val{:integrator}},
                           tangent_prob::DEProblem,
                           tangent_state = (@view tangent_prob.u0[:, 2:end]);
                           phase_state = (@view tangent_prob.u0[:, 1]),
                           kwargs...)
    u0 = augmented_vector(phase_state, tangent_state)
    new_prob = with_u0(tangent_prob, u0)
    integrator = get_integrator(new_prob; kwargs...)
    keepgoing!(integrator)
    return integrator
end

tangent_propagate(::Type{Val{:integrator}}, le_solver::LESolver,
                  tangent_state = le_solver.tangent_state;
                  phase_state = le_solver.phase_state,
                  kwargs...) =
    tangent_propagate(Val{:integrator}, de_prob(le_solver), tangent_state;
                      phase_state = phase_state,
                      kwargs...)
