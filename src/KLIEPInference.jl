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
  spKLIEP, spKLIEP!,
  Hinv_row,

  # inference
  BootstrapEstimates,
  boot_SparKLIE,
  simulCI, simulCIstudentized,

  # other
  Ψising, unpack, pack,
  KLIEP_Hessian,
  SparKLIE_stderr

include("utils.jl")
include("bootstrap.jl")
include("sampler.jl")
include("solver.jl")
include("debias.jl")
include("stderr.jl")

end # module
