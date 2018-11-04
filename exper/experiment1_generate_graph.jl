using KLIEPInference
using Distributions
using LinearAlgebra, Random
using JLD

p = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])

θx = chain(p, 5, 0.2, 0.4, sgn)
θy = copy(θx)
Δ = θx - θy

@save "params_exp1_$(p)_$(sgn).jld" p sgn θx θy
