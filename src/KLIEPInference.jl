module KLIEPInference


using Distributions, LinearAlgebra, Statistics
import JuMP, MathOptInterface, MathOptInterfaceMosek, SCS
const MOI = MathOptInterface
using StatsBase
using ProximalBase, CoordinateDescent

export
  IsingSampler,
  chain,

  # KLIEP solvers
  CD_KLIEP,
  SCS_KLIEP,
  Mosek_KLIEP,
  KLIEP, KLIEP!,
  spKLIEP, spKLIEP!,
  Hinv_row,

  # inference
  boot_KLIEP,
  simulCI,
  simulCIstudentized,

  # utils
  Ψising,
  KLIEP_Hessian

include("utils.jl")
include("sampler.jl")
include("solver.jl")
include("bootstrap.jl")





end # module
