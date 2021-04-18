using KLIEPInference
using LinearAlgebra, Random
using Distributions, JLD

Random.seed!(543)

m = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])

γx = chain(m; lenC=5, lb=0.2, ub=0.4, sgn=sgn)
γy = copy(γx)

@save "graphs/params_exp3_$(m)_$(sgn).jld" m sgn γx γy