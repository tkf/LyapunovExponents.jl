"""
A type to hold an example dynamical system and its known Lyapunov exponents.

Here is an example code for constructing an example dynamical system,
calculate its LEs and plot them:
```julia
using LyapunovExponents
using Plots
demo = solve!(LEDemo(LyapunovExponents.lorenz_63()))
plot(demo)
```
"""
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
Return a [`LEExample`](@ref) for the Lorenz system.

* <https://en.wikipedia.org/wiki/Lorenz_system>
* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* E. N. Lorenz, J. Atmos. Sci. 20, 130-141 (1963)
"""
function lorenz_63(;
        u0=[0.1, 0.1, 0.1],
        tspan=(0.0, 1.0),
        num_attr=4000,
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
Return a [`LEExample`](@ref) for the simplest piecewise linear
dissipative chaotic flow.

* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
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
Return a [`LEExample`](@ref) for the van der Pol oscillator with
periodic forcing.

`.known_exponents` are extracted from Figure 6 of Geist, Parlitz &
Lauterborn (1990).

* <http://scholarpedia.org/article/Van_der_Pol_oscillator>
* <https://en.wikipedia.org/wiki/Van_der_Pol_oscillator>
* van der Pol and van der Mark. “Frequency Demultiplication.”
  Nature 120, no. 3019 (September 1927): 363.
  https://doi.org/10.1038/120363a0.
* Parlitz, Ulrich, and Werner Lauterborn. “Period-Doubling Cascades and
  Devil’s Staircases of the Driven van Der Pol Oscillator.”
  Physical Review A 36, no. 3 (August 1, 1987): 1428–34.
  https://doi.org/10.1103/PhysRevA.36.1428.
  (Figure 10a)
* Geist, K., Parlitz, U., & Lauterborn, W. (1990).
  Comparison of Different Methods for Computing Lyapunov Exponents.
  Progress of Theoretical Physics, 83, 875–893.
  https://doi.org/10.1143/PTP.83.875
  (Figure 6)
"""
function van_der_pol(;
        u0=[0.1, 0.1],
        t0=0.0, chunk_periods=1,
        ω=2.466, tspan=(t0, t0 + chunk_periods * 2 * π / ω),
        num_attr=100,
        atol=0, rtol=1e-1)
    # Note that with larger num_attr (e.g., 10000), the last Lyapunov
    # exponents negatively overshoots what Geist, Parlitz & Lauterborn
    # (1990) reported.  num_attr=100 is required for test to pass.
    @inline function phase_dynamics!(t, u, du)
        du[1] = u[2]
        du[2] = -u[1] - 5.0 * (u[1]^2 - 1) * u[2] + 5.0 * cos(ω * t)
    end
    ContinuousExample(
        "van der Pol & van der Mark (1927)",
        phase_dynamics!,
        u0, tspan, num_attr,
        [0.085, -6.7],   # known_exponents
        atol, rtol,
    )
end

σ(x) = 1 / (1 + exp(-x))

type ContinuousRNN
    w
    θ
    τ
end

function (param::ContinuousRNN)(t, u, du)
    w = param.w
    θ = param.θ
    τ = param.τ
    du .= (- u .+ w * σ.(u .+ θ)) ./ τ
end

"""
Return a [`LEExample`](@ref) for a low-dimensional chaotic
continuous-time recurrent neural networks by Beer (1995).

* Beer, R. D. (1995). On the dynamics of small continuous-time recurrent
  neural networks. Adapt. Behav., 3(4), 469–509.
  https://doi.org/10.1177/105971239500300405
  (Figure 9D)
"""
function beer_95(;
        u0=[0.1, 0.1, 0.1],
        tspan=(0, 10.0),
        num_attr=4000,
        atol=0, rtol=1e-1)
    phase_dynamics! = ContinuousRNN(
        # w
        [  5.422  -0.018  2.75
          -0.24    4.59   1.21
           0.535  -2.25   3.885 ],
        # θ
        [-4.108, -2.787, -1.114],
        # τ
        [1.0, 2.5, 1.0],
    )
    ContinuousExample(
        "Beer (1995)",
        phase_dynamics!,
        u0, tspan, num_attr,
        [0.010, 0],   # known_exponents
        atol, rtol,
    )
end

"""
Return a [`LEExample`](@ref) for the Hénon map.

* M. Hénon, Commun. Math. Phys. Phys. 50, 69-77 (1976)
* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* <https://en.wikipedia.org/wiki/H%C3%A9non_map>
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
Return a [`LEExample`](@ref) for the Chirikov standard map.

* B. V. Chirikov, Physics Reports 52, 263-379 (1979)
* <http://sprott.physics.wisc.edu/chaos/comchaos.htm>
* <https://en.wikipedia.org/wiki/Standard_map>
* <http://www.scholarpedia.org/article/Chirikov_standard_map>
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

* <https://en.wikipedia.org/wiki/Baker%27s_map>
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
Return a [`LEExample`](@ref) for the Arnold's cat map

* <https://en.wikipedia.org/wiki/Arnold%27s_cat_map>
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
    van_der_pol,
    beer_95,
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

"""
    LEDemo(example::LEExample; <keyword arguments>)

Create a `LEDemo` holding an `example` and an appropriate `LEProblem`
created from the `example`.
"""
function LEDemo(example::LEExample; kwargs...)
    LEDemo(example, le_problem(example; kwargs...))
end

"""
    solve!(demo::LEDemo; progress=-1, <keyword arguments>)

Initialize `demo.solver` from `demo.prob` and run
`solve!(demo.solver)` to calculate the Lyapunov exponents.
"""
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
