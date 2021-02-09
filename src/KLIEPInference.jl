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
  CD_KLIEP,
  KLIEP, KLIEP!,
  spKLIEP, spKLIEP!, spKLIEP_refit!,
  Hinv_row,

  # inference
  BootstrapEstimates,
  boot_SparKLIE1,
  boot_KLIEP,
  boot_spKLIEP, boot_spKLIEPfull,
  boot_oracleKLIEP,
  boot_gaussKLIEP,
  simulCI,
  simulCIstudentized,

  # other
  Ψising, unpack, pack,
  KLIEP_Hessian,
  KLIEP_debias1, KLIEP_debias2,
  KLIEP_stderr

include("utils.jl")
include("bootstrap.jl")
include("sampler.jl")
include("solver.jl")
include("debias.jl")
include("stderr.jl")

end # module
