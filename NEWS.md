# Unreleased

Use CoordinateDescent 0.8's Julia-style API (`coordinate_descent!`,
`num_coordinates`, `descend_coordinate!`, `warm_start`) and ProximalBase 0.5's
`NormL1`, `ZeroFunction` and `prox_coordinate!`. The KLIEP fitting functions
continue to return the supplied iterate by extracting `CDResults.minimizer`.
Julia 1.10 or later is required. The root and test environments now resolve
from Project.toml rather than the committed Julia 1.5 manifests.

The source overrides use the dependencies' masters on GitHub until
CoordinateDescent 0.8 and its solver dependencies are registered.
