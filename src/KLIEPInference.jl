module KLIEPInference

using Distributions, LinearAlgebra, Statistics, SparseArrays, Random
using StatsBase
using ProximalBase, CoordinateDescent

export
    IsingSampler,
    chain,
    kary_tree,
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
    boot_quantile,
    boot_quantile_studentized,

    # other
    Ψising,
    unpack,
    pack,
    KLIEP_Hessian,
    stderr_SparKLIE

include("utils.jl")
include("bootstrap.jl")
include("sampler.jl")
include("solver.jl")
include("debias.jl")
include("stderr.jl")

end # module
