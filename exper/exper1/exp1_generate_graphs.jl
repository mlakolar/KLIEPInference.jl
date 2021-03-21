using KLIEPInference
using Random
using Distributions
using JLD

Random.seed!(325)

## chain (1)
idx0 = KLIEPInference.trimap(5, 6)
idx1 = KLIEPInference.trimap(4, 5)
idx2 = KLIEPInference.trimap(6, 7)
idx3 = KLIEPInference.trimap(4, 6)
idx4 = KLIEPInference.trimap(5, 7)

for m in [25, 50]
    γx = chain(m; lenC=m, sgn=0)
    γy = deepcopy(γx)
    γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
    γy[idx1] *= -1.0
    γy[idx2] *= -1.0
    γy[idx3] += (2. * rand(Bernoulli()) - 1.) * 0.1
    γy[idx4] += (2. * rand(Bernoulli()) - 1.) * 0.1
    @save "graphs/params_exp1_chain1_$(m).jld" γx γy
end

## chain (2)
idx0 = KLIEPInference.trimap(5, 6)
idx1 = KLIEPInference.trimap(3, 4)
idx2 = KLIEPInference.trimap(7, 8)
idx3 = KLIEPInference.trimap(4, 5)
idx4 = KLIEPInference.trimap(4, 6)

for m in [25, 50]
    γx = chain(m; lenC=m, sgn=0)
    γy = deepcopy(γx)
    γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
    γy[idx1] *= -1.0
    γy[idx2] *= -1.0
    γy[idx3] = abs(γy[idx3]) > 0.2 ? γy[idx3] - sign(γy[idx3])*0.1 : γy[idx3] + sign(γy[idx3])*0.1
    γy[idx4] += (2. * rand(Bernoulli()) - 1.) * 0.1
    @save "graphs/params_exp1_chain2_$(m).jld" γx γy
end

Random.seed!(324)

## tree (1)
idx0 = KLIEPInference.trimap(1, 3)
idx1 = KLIEPInference.trimap(1, 2)
idx2 = KLIEPInference.trimap(3, 8)
idx3 = KLIEPInference.trimap(1, 9)
idx4 = KLIEPInference.trimap(3, 4)

for m in [25, 50]
    γx = kary_tree(m, sgn=0)
    γy = deepcopy(γx)
    γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
    γy[idx1] *= -1.0
    γy[idx2] *= -1.0
    γy[idx3] += (2. * rand(Bernoulli()) - 1.) * 0.1
    γy[idx4] = abs(γy[idx4]) > 0.2 ? γy[idx4] - sign(γy[idx4])*0.1 : γy[idx4] + sign(γy[idx4])*0.1
    @save "graphs/params_exp1_tree1_$(m).jld" γx γy
end

## tree (2)
idx0 = KLIEPInference.trimap(1, 3)
idx1 = KLIEPInference.trimap(2, 5)
idx2 = KLIEPInference.trimap(8, 23)
idx3 = KLIEPInference.trimap(2, 3)
idx4 = KLIEPInference.trimap(3, 4)

for m in [25, 50]
    γx = kary_tree(m, sgn=0)
    γy = deepcopy(γx)
    γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
    γy[idx1] *= -1.0
    γy[idx2] *= -1.0
    γy[idx3] += (2. * rand(Bernoulli()) - 1.) * 0.1
    γy[idx4] += (2. * rand(Bernoulli()) - 1.) * 0.1
    @save "graphs/params_exp1_tree2_$(m).jld" γx γy
end
