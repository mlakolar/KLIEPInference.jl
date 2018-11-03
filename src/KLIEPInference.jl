module KLIEPInference


using Distributions, LinearAlgebra, Statistics
import JuMP, MathOptInterface, MathOptInterfaceMosek, SCS
const MOI = MathOptInterface
using StatsBase
using ProximalBase, CoordinateDescent

export
  IsingSampler,

  # KLIEP solvers
  CD_KLIEP,
  SCS_KLIEP,
  Mosek_KLIEP,
  KLIEP, KLIEP!,
  spKLIEP,
  Hinv_row,

  # inference
  boot_KLIEP,
  simulCI,
  simulCIstudentized,

  # utils
  Ψising,
  KLIEP_Hessian


include("sampler.jl")
include("solver.jl")
include("utils.jl")
include("bootstrap.jl")





end # module
