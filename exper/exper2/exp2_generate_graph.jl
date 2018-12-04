using KLIEPInference
using Distributions
using LinearAlgebra, Random
using JLD

Random.seed!(543)

p = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
numChanges = parse(Int,ARGS[3])
lbInd = parse(Int,ARGS[4])

lbArr = 0.0:.05:0.5
lb = lbArr[lbInd]

θx = chain(p, 5, 0.2, 0.4, sgn)
θy = copy(θx)


ind_change = sample(1:length(θx), numChanges; replace=false)
du = Uniform(lb, lb+0.01)
for i in ind_change
  θy[i] -= rand(du)
end

Δ = θx - θy
@show ind_nz = findall(!iszero, Δ)
@show Δ[ind_nz]

@save "/home/byolkim/graphs/params_exp2_$(p)_$(sgn)_$(numChanges)_$(lbInd).jld" p sgn θx θy
