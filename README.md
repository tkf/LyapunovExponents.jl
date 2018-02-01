# LyapunovExponents.jl --- A hackable Lyapunov exponents calculator

[![Build Status](https://travis-ci.org/tkf/LyapunovExponents.jl.svg?branch=master)](https://travis-ci.org/tkf/LyapunovExponents.jl)
[![Coverage Status](https://coveralls.io/repos/tkf/LyapunovExponents.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/tkf/LyapunovExponents.jl?branch=master)
[![codecov.io](http://codecov.io/github/tkf/LyapunovExponents.jl/coverage.svg?branch=master)](http://codecov.io/github/tkf/LyapunovExponents.jl?branch=master)
[![](https://img.shields.io/badge/docs-stable-red.svg)](https://tkf.github.io/LyapunovExponents.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://tkf.github.io/LyapunovExponents.jl/latest)

The aim of LyapunovExponents.jl is to provide an efficient research
platform for computations related to Lyapunov exponents.  This is
(planned to be) achieved by exposing low-level APIs to Lyapunov
exponents calculation.

At the moment, LyapunovExponents.jl is still at the very early stage
of development and nowhere close to providing a stable API.


## Features

Implemented:

* Lyapunov exponents calculation based on QR decomposition.
* Maximum Lyapunov exponent calculation.
* Covariant Lyapunov vectors calculation.
  (Not well tested!)
* Tangent space evolution based on the automatic differentiation tool
  [ForwardDiff.jl].
* Testing utilities for tangent space evolution (Jacobian calculation)
  provided by users.
* Continuous dynamical systems support based on [DifferentialEquations.jl].
  This means that [the rich set of ODE solvers](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html)
  can be used.
* [Various examples of continuous and discrete dynamical systems.](https://tkf.github.io/LyapunovExponents.jl/latest/examples/)
* [OnlineStats.jl] support: on-the-fly calculation of Lyapunov
  exponents, their variance, covariance, and any other statics it
  supports.

Wanted list:

* Delay differential equations.
* Partial differential equations.
* Stochastic dynamical systems.
* Poincaré map.

[DifferentialEquations.jl]: http://juliadiffeq.org
[ForwardDiff.jl]: http://www.juliadiff.org/ForwardDiff.jl
[ChaosTools.jl]: https://juliadynamics.github.io/DynamicalSystems.jl/latest/chaos/overview/
[DynamicalSystems.jl]: https://juliadynamics.github.io/DynamicalSystems.jl/latest/
[OnlineStats.jl]: https://github.com/joshday/OnlineStats.jl
[eom]: https://github.com/termoshtt/eom


## Related works

* [ChaosTools.jl] from the [DynamicalSystems.jl] ecosystem is another
  Julia library which provides easy-to-use, clearly written,
  well-tested, and well-documented Lyapunov exponents calculation.
* [eom] is a Rust library which provides Lyapunov exponents and
  Covariant Lyapunov vectors calculation, on top of configurable
  ODE/PDE solvers.
* (I'm sure there are more...)


## License

The LyapunovExponents.jl package is licensed under the MIT "Expat" License.
See [LICENSE.md]() file.
