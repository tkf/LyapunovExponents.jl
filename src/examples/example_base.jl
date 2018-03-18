module ExampleBase

import DifferentialEquations: solve, solve!
using Parameters: @with_kw

using ...LyapunovExponents: LEProblem, ContinuousLEProblem, DiscreteLEProblem,
    init
import ...LyapunovExponents: dimension, report

"""
A type to hold an example dynamical system and its known Lyapunov exponents.
"""
@with_kw struct LEExample{ProblemType}
    name
    phase_dynamics
    u0
    t_renorm
    param
    tangent_dynamics
    t_attr
    known_exponents
    atol
    rtol
    terminator_options
    de_options = []
end

LEExample{ProblemType}(
            name,
            phase_dynamics,
            u0,
            t_renorm,
            param,
            tangent_dynamics,
            t_attr,
            known_exponents,
            atol,
            rtol,
            terminator_options;
            kwargs...
        ) where {ProblemType} =
    LEExample{ProblemType}(;
            name = name,
            phase_dynamics = phase_dynamics,
            u0 = u0,
            t_renorm = t_renorm,
            param = param,
            tangent_dynamics = tangent_dynamics,
            t_attr = t_attr,
            known_exponents = known_exponents,
            atol = atol,
            rtol = rtol,
            terminator_options = terminator_options,
            kwargs...)

const ContinuousExample = LEExample{ContinuousLEProblem}
const DiscreteExample = LEExample{DiscreteLEProblem}

function LEProblem(example::LEExample{P}; kwargs...) where {P <: LEProblem}
    t_tran = example.t_attr * 0.1
    if P <: DiscreteLEProblem
        t_tran = ceil(typeof(example.t_attr), t_tran)
    end
    P(example.phase_dynamics,
      example.u0,
      example.param;
      t_renorm = example.t_renorm,
      t_attr = example.t_attr,
      t_tran = t_tran,
      tangent_dynamics = example.tangent_dynamics,
      kwargs...)
end

dimension(example::LEExample) = length(example.u0)

function solve(example::LEExample;
               dim_lyap=dimension(example), kwargs...)
    solve(LEProblem(example; dim_lyap=dim_lyap);
          terminator_options = example.terminator_options,
          de_options = example.de_options,
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
    demo.solver = solve(demo.prob;
                        progress = progress,
                        record = record,
                        terminator_options = demo.example.terminator_options,
                        de_options = demo.example.de_options,
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

function report(io::IO, demo::LEDemo)
    print_with_color(:blue, io, "Lyapunov Exponents Demo")
    println(io, ": ", demo.example.name)

    if isdefined(demo, :solver)
        report(io, demo.solver)
    else
        print_with_color(:red, io, "[solver not initialized]")
        println(io)
    end

    if ! isempty(demo.example.known_exponents)
        print_with_color(:yellow, io, "Known LEs")
        print(io, ": ")
        show(IOContext(io, :limit => true), demo.example.known_exponents)
        println(io)
    end
end

end
