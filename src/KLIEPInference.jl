module KLIEPInference


using Distributions, LinearAlgebra, Statistics, SparseArrays
import JuMP, MathOptInterface, MathOptInterfaceMosek, SCS
const MOI = MathOptInterface
using StatsBase
using ProximalBase, CoordinateDescent

export
  IsingSampler,
  chain,
  removeEdges!,
  addEdges!,

  # KLIEP solvers
  CD_KLIEP,
  SCS_KLIEP,
  Mosek_KLIEP,
  KLIEP, KLIEP!,
  spKLIEP, spKLIEP!, spKLIEP_refit!,
  Hinv_row,

  # inference
  boot_KLIEP,
  boot_spKLIEP,
  simulCI,
  simulCIstudentized,

  # utils
  Ψising,
  KLIEP_Hessian,
  unpack

include("utils.jl")
include("sampler.jl")
include("solver.jl")
include("bootstrap.jl")





end # module
