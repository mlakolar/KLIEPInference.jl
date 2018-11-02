module KLIEPInference


using Distributions, LinearAlgebra
import JuMP, MathOptInterface, MathOptInterfaceMosek, SCS
const MOI = MathOptInterface
using ProximalBase, CoordinateDescent

export
  IsingSampler,

  # KLIEP solvers
  CD_KLIEP,
  SCS_KLIEP,
  Mosek_KLIEP,
  KLIEP,
  spKLIEP,

  #
  Ψising


include("sampler.jl")
include("solver.jl")
include("utils.jl")





end # module
