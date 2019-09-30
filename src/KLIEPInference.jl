module KLIEPInference

using Distributions, LinearAlgebra, Statistics, SparseArrays, Random
using StatsBase
using ProximalBase, CoordinateDescent

export
  IsingSampler,
  chain,
  removeEdges!,
  addEdges!,

  # KLIEP solvers
  KLIEPSolver,
  CD_KLIEP, CD_SymKLIEP,
  KLIEP, KLIEP!, SymKLIEP, SymKLIEP!,
  spKLIEP, spKLIEP!, spKLIEP_refit!,
  spSymKLIEP, spSymKLIEP!, spSymKLIEP_refit!,
  Hinv_row,

  # inference
  BootstrapEstimates,
  boot_KLIEP,
  boot_spKLIEP, boot_spKLIEPfull,
  boot_oracleKLIEP,
  boot_gaussKLIEP,
  simulCI,
  simulCIstudentized,

  # utils
  Ψising,
  KLIEP_Hessian, SymKLIEP_Hessian,
  KLIEP_debias, SymKLIEP_debias,
  KLIEP_var, SymKLIEP_var,
  unpack

include("utils.jl")
include("utils_single_edge.jl")
include("utils_SymKLIEPLoss.jl")
include("sampler.jl")
include("solver.jl")
include("solver_SymKLIEP.jl")
include("bootstrap.jl")



end # module
