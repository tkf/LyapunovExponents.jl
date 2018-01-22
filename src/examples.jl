struct LEExample{ProblemType}
    name
    phase_dynamics!
    u0
    tspan
    num_attr
    known_exponents
    atol
    rtol
end

const ContinuousExample = LEExample{ContinuousLEProblem}
const DiscreteExample = LEExample{DiscreteLEProblem}

"""
Lorenz system.

* https://en.wikipedia.org/wiki/Lorenz_system
* http://sprott.physics.wisc.edu/chaos/comchaos.htm
* E. N. Lorenz, J. Atmos. Sci. 20, 130-141 (1963)
"""
function lorenz_63(;
        u0=[0.1, 0.1, 0.1],
        tspan=(0.0, 1.0),
        num_attr=1000,
        atol=0, rtol=1e-2)
    @inline function phase_dynamics!(t, u, du)
        du[1] = 10.0(u[2]-u[1])
        du[2] = u[1]*(28.0-u[3]) - u[2]
        du[3] = u[1]*u[2] - (8/3)*u[3]
    end
    ContinuousExample(
        "Lorenz (1963)",
        phase_dynamics!,
        u0, tspan, num_attr,
        [0.9056, 0, -14.5723],  # known_exponents
        atol, rtol,
    )
end

"""
Simplest piecewise linear dissipative chaotic flow.

* http://sprott.physics.wisc.edu/chaos/comchaos.htm
* S. J. Linz and J. C. Sprott, Phys. Lett. A 259, 240-245 (1999)
"""
function linz_sprott_99(;
        u0=[0.1, 0.1, 0.1],
        tspan=(0.0, 1.0),
        num_attr=10000,
        atol=0, rtol=1e-2)
    @inline function phase_dynamics!(t, u, du)
        du[1] = u[2]
        du[2] = u[3]
        du[3] = -0.6 * u[3] - u[2] - (u[1] > 0 ? u[1] : -u[1]) + 1
    end
    ContinuousExample(
        "Linz & Sprott (1999) Piecewise linear flow",
        phase_dynamics!,
        u0, tspan, num_attr,
        [0.0362, 0, -0.6362],   # known_exponents
        atol, rtol,
    )
end

"""
Hénon map.

* M. Hénon, Commun. Math. Phys. Phys. 50, 69-77 (1976)
* http://sprott.physics.wisc.edu/chaos/comchaos.htm
* https://en.wikipedia.org/wiki/H%C3%A9non_map
"""
function henon_map(;
        u0=[0.1, 0.1],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=1e-2)
    @inline function phase_dynamics!(t, u, u_next)
        u_next[1] = 1 + u[2] - 1.4 * u[1]^2
        u_next[2] = 0.3 * u[1]
    end
    DiscreteExample(
        "Hénon map",
        phase_dynamics!,
        u0, tspan, num_attr,
        [0.41922, -1.62319],   # known_exponents
        atol, rtol,
    )
end

"""
Chirikov standard map.

* B. V. Chirikov, Physics Reports 52, 263-379 (1979)
* http://sprott.physics.wisc.edu/chaos/comchaos.htm
* https://en.wikipedia.org/wiki/Standard_map
* http://www.scholarpedia.org/article/Chirikov_standard_map
"""
function standard_map(;
        u0=[2.68156, 2.31167],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=0.2)
    # TODO: Improve the accuracy. Check the paper.  It looks like
    # `num_attr=1000000` is required to see some kind of convergence.
    @inline function phase_dynamics!(t, u, u_next)
        u_next[2] = (u[2] + sin(u[1])) % 2π
        u_next[1] = (u[1] + u_next[2]) % 2π
    end
    DiscreteExample(
        "Chirikov standard map",
        phase_dynamics!,
        u0, tspan, num_attr,
        [0.10497, -0.10497],   # known_exponents
        atol, rtol,
    )
end

"""
Baker's map

* https://en.wikipedia.org/wiki/Baker%27s_map
"""
function bakers_map(;
        u0=[0.6, 0.4],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=1e-7)
    @inline function phase_dynamics!(t, u, u_next)
        if u[1] < 0.5
            u_next[1] = 2 * u[1]
            u_next[2] = u[2] / 2
        else
            u_next[1] = 2 - 2 * u[1]
            u_next[2] = 1 - u[2] / 2
        end
    end
    DiscreteExample(
        "Baker's map",
        phase_dynamics!,
        u0, tspan, num_attr,
        log.([2.0, 0.5]),       # known_exponents
        atol, rtol,
    )
end

"""
Arnold's cat map

* https://en.wikipedia.org/wiki/Arnold%27s_cat_map
"""
function arnold_cat_map(;
        u0=[0.1, 0.1],
        tspan=(0, 10),
        num_attr=10000,
        atol=0, rtol=1e-4)
    M = [2 1; 1 1]
    @inline function phase_dynamics!(t, u, u_next)
        A_mul_B!(u_next, M, u)
        u_next .= mod.(u_next, 1)
    end
    DiscreteExample(
        "Arnold's cat map",
        phase_dynamics!,
        u0, tspan, num_attr,
        log.(sort!(eigvals(M), rev=true)), # known_exponents
        atol, rtol,
    )
end

const EXAMPLES = [
    lorenz_63,
    linz_sprott_99,
    henon_map,
    standard_map,
    bakers_map,
    arnold_cat_map,
]

function le_problem(example::ContinuousExample; kwargs...)
    ContinuousLEProblem(
        example.phase_dynamics!,
        example.u0,
        example.tspan;
        kwargs...)
end

function le_problem(example::DiscreteExample; kwargs...)
    DiscreteLEProblem(
        example.phase_dynamics!,
        example.u0,
        example.tspan;
        kwargs...)
end

dimension(example::LEExample) = length(example.u0)

function solve(example::LEExample;
               dim_lyap=dimension(example), kwargs...)
    solve(le_problem(example; dim_lyap=dim_lyap), example.num_attr;
          kwargs...)
end

mutable struct LEDemo
    example
    prob
    solver

    LEDemo(example, prob) = new(example, prob)
end

function LEDemo(example::LEExample; kwargs...)
    LEDemo(example, le_problem(example; kwargs...))
end

function solve!(demo::LEDemo; progress=-1, kwargs...)
    demo.solver = LERecordingSolver(init(demo.prob; progress=progress),
                                    demo.example.num_attr)
    solve!(demo.solver; progress=progress, kwargs...)
end

function Base.show(io::IO, demo::LEDemo)
    print(io, "Demo: ", demo.example.name)
    if isdefined(demo, :solver)
        print(io, ", ", demo.solver)
    else
        print(io, " [solver not initialized]")
    end
end
