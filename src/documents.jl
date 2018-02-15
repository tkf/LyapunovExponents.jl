import DiffEqBase: init
"""
    init(prob::LEProblem; <keyword arguments>) :: LESolver
    init(prob::CLVProblem; <keyword arguments>) :: CLVSolver

These are simply the aliases of `LESolver(prob; <keyword arguments>)`
and `CLVSolver(prob; <keyword arguments>)`.  See [`LESolver`](@ref)
and [`CLVSolver`](@ref) for supported keyword arguments.
"""
function init end

import DiffEqBase: solve
"""
    solve(prob::LEProblem; <keyword arguments>) :: LESolution
    solve(prob::CLVProblem; <keyword arguments>) :: CLVSolution

Equivalent to:

```julia
solver = init(prob; <keyword arguments except progress>)
solve!(solver; progress=progress)
solver.sol
```
"""
function solve end
using DiffEqBase: solve

import DiffEqBase: solve!
"""
    solve!(solver::LESolver; <keyword arguments>) :: LESolver
    solve!(solver::CLVSolver; <keyword arguments>) :: CLVSolver

Solve pre-initialized problem.
"""
function solve! end
using DiffEqBase: solve!
