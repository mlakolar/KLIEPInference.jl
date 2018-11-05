using KLIEPInference
using Distributions
using LinearAlgebra, Random
using JLD

p = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
rem = parse(Int,ARGS[3])

θx = chain(p, 5, 0.2, 0.4, sgn)
θy = copy(θx)

removeEdges!(θy, rem)
Δ = θx - θy
@show ind_nz = findall(!iszero, Δ)
@show Δ[ind_nz]

@save "params_exp2_$(p)_$(sgn).jld" p sgn θx θy
