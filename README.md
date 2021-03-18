# KLIEPInference.jl

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

You will also need to install the following packages.

- JLD, Distributions, StatsBase, Pyplot, IJulia

`using Pkg; Pkg.add.(["JLD", "Distributions", "StatsBase", "Pyplot", "IJulia"])`