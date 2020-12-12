module KLIEPInference

using Distributions, LinearAlgebra, Statistics, SparseArrays, Random
using StatsBase
using ProximalBase, CoordinateDescent
using Printf

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
  Hinv_row, Hinv_row_refit!,

  # inference
  BootstrapEstimates,
  boot_KLIEP,
  boot_spKLIEP, boot_spKLIEPfull,
  boot_oracleKLIEP,
  boot_gaussKLIEP,
  simulCI,
  simulCIstudentized,

  # other
  Ψising, unpack, pack,
  KLIEP_Hessian,
  debias_KLIEP, stderr_KLIEP, tuneλ

include("utils.jl")
include("bootstrap.jl")
include("sampler.jl")
include("solver.jl")
include("debias.jl")
include("stderr.jl")
include("tune.jl")

end # module
