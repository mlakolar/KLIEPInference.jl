using KLIEPInference
using LinearAlgebra, Random
using Distributions, JLD

Random.seed!(543)

m = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
numChanges = parse(Int,ARGS[3])
lbInd = parse(Int,ARGS[4])

lbArr = 0.0:.05:0.5
lb = lbArr[lbInd]

γx = chain(m, 5, 0.2, 0.4, sgn)
γy = copy(γx)

ind_change = sample(1:length(γx), numChanges; replace=false)
du = Uniform(lb, lb+0.01)
for i in ind_change
    γy[i] -= rand(du)
end

θ = γx - γy
@show ind_nz = findall(!iszero, θ)
@show θ[ind_nz]

@save "params_exp4_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld" m sgn γx γy