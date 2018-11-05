module KLIEPInference

using Distributions, LinearAlgebra, Statistics, SparseArrays
using StatsBase
using ProximalBase, CoordinateDescent

export
  IsingSampler,
  chain,
  removeEdges!,
  addEdges!,

  # KLIEP solvers
  KLIEPSolver,
  CD_KLIEP,
  KLIEP, KLIEP!,
  spKLIEP, spKLIEP!, spKLIEP_refit!,
  Hinv_row,

  # inference
  boot_KLIEP,
  boot_spKLIEP,
  boot_oracleKLIEP,
  boot_gaussKLIEP,
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
