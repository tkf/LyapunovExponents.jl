module ExampleBase

import DifferentialEquations: solve, solve!

using ...LyapunovExponents: LEProblem, ContinuousLEProblem, DiscreteLEProblem
import ...LyapunovExponents: dimension

"""
A type to hold an example dynamical system and its known Lyapunov exponents.
"""
struct LEExample{ProblemType}
    name
    phase_dynamics!
    u0
    tspan
    num_attr
    tangent_dynamics!
    known_exponents
    atol
    rtol
end

const ContinuousExample = LEExample{ContinuousLEProblem}
const DiscreteExample = LEExample{DiscreteLEProblem}

function LEProblem(example::LEExample{P}; kwargs...) where {P <: LEProblem}
    P(example.phase_dynamics!,
      example.u0,
      example.tspan;
      tangent_dynamics! = example.tangent_dynamics!,
      kwargs...)
end

dimension(example::LEExample) = length(example.u0)

function solve(example::LEExample;
               dim_lyap=dimension(example), kwargs...)
    solve(LEProblem(example; dim_lyap=dim_lyap), example.num_attr;
          kwargs...)
end

mutable struct LEDemo
    example
    prob
    solver

    LEDemo(example, prob) = new(example, prob)
end

"""
    LEDemo(example::LEExample; <keyword arguments>)

Here is an example code for constructing an example dynamical system,
calculate its LEs and plot them:
```julia
using LyapunovExponents
using Plots
demo = solve!(LyapunovExponents.lorenz_63())
plot(demo)
```

Create a `LEDemo` holding an `example` and an appropriate `LEProblem`
created from the `example`.
"""
function LEDemo(example::LEExample; kwargs...)
    LEDemo(example, LEProblem(example; kwargs...))
end

"""
    solve!(demo::LEDemo; progress=-1, <keyword arguments>)

Initialize `demo.solver` from `demo.prob` and run
`solve!(demo.solver)` to calculate the Lyapunov exponents.
"""
function solve!(demo::LEDemo; progress = -1, record = true, kwargs...)
    demo.solver = solve(demo.prob, demo.example.num_attr;
                        progress = progress,
                        record = record,
                        kwargs...)
    return demo
end

function Base.show(io::IO, demo::LEDemo)
    print(io, "Demo: ", demo.example.name)
    if isdefined(demo, :solver)
        print(io, ", ", demo.solver)
    else
        print(io, " [solver not initialized]")
    end
end

end
