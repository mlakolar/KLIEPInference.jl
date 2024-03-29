# KLIEPInference.jl

[![Testing](https://github.com/mlakolar/KLIEPInference.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mlakolar/KLIEPInference.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mlakolar/KLIEPInference.jl/branch/master/graph/badge.svg?token=XVJH400GOY)](https://codecov.io/gh/mlakolar/KLIEPInference.jl)

Code for the paper
_Two-sample inference for high-dimensional Markov networks_
by _Byol Kim_, _Song Liu_, and  _Mladen Kolar_. [arXiv: 1905.00466](https://arxiv.org/abs/1905.00466)

The code has been tested with Julia 1.5.4.

Install Julia for your platform: [https://julialang.org/downloads/](https://julialang.org/downloads/).

### Installing the package

You can obtain KLIEPInference using Julia's Pkg REPL-mode (hitting `]` as the first character of the command prompt):
```
(@v1.5) pkg> add https://github.com/mlakolar/KLIEPInference.jl
```

You will also need to install the following software.

- JLD, Distributions, ProximalBase, CoordinateDescent, StatsBase, PyPlot, IJulia, Revise.

`using Pkg; Pkg.add.(["JLD", "Distributions", "ProximalBase", "CoordinateDescent", "StatsBase", "PyPlot", "IJulia", "Revise"])`

- GNU Parallel.

You can install this package using `sudo apt install parallel`.

- Jupyter Notebook.

To install Jupyter Notebook, follow instructions [here](https://jupyter.org/install).

### Reproducing experimental results

Please read "README" at each experiment folder. 
