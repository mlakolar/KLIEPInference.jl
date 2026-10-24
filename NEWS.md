# Unreleased

Use CoordinateDescent 0.8's Julia-style API (`coordinate_descent!`,
`num_coordinates`, `descend_coordinate!`, `warm_start`) and ProximalBase 0.5's
`NormL1`, `ZeroFunction` and `prox_coordinate!`. The KLIEP fitting functions
continue to return the supplied iterate by extracting `CDResults.minimizer`.
Julia 1.10 or later is required. The root and test environments now resolve
from Project.toml rather than the committed Julia 1.5 manifests.

The CoordinateDescent source override temporarily points to the branch for
mlakolar/CoordinateDescent.jl#40; merge that change before this update and
switch the override to master or a registered 0.8 release.
