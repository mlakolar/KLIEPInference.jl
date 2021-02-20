using KLIEPInference
using Distributions
using LinearAlgebra, Random
using JLD

Random.seed!(543)

m = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])

θx = chain(m, 5, 0.2, 0.4, sgn)
θy = copy(θx)

@save "params_exp3_$(m)_$(sgn).jld" m sgn θx θy
