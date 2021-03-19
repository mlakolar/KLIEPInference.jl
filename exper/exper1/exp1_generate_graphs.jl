using KLIEPInference
using Random
using JLD

Random.seed!(325)

## chain (1)
idx0 = KLIEPInference.trimap(5, 6)
idx1 = KLIEPInference.trimap(4, 5)
idx2 = KLIEPInference.trimap(6, 7)
idx3 = KLIEPInference.trimap(4, 6)
idx4 = KLIEPInference.trimap(5, 7)

# m = 25
γx = chain(25, 25, 0.1, 0.3, 0)
γy = deepcopy(γx)
γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
γy[idx1] *= -1.0
γy[idx2] *= -1.0
γy[idx3] -= 0.1
γy[idx4] -= 0.1
@save "graphs/params_exp1_chain1_25.jld" γx γy

# m = 50
γx = chain(50, 50, 0.1, 0.3, 0)
γy = deepcopy(γx)
γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
γy[idx1] *= -1.0
γy[idx2] *= -1.0
γy[idx3] -= 0.1
γy[idx4] -= 0.1
@save "graphs/params_exp1_chain1_50.jld" γx γy

## chain (2)
idx0 = KLIEPInference.trimap(5, 6)
idx1 = KLIEPInference.trimap(3, 4)
idx2 = KLIEPInference.trimap(7, 8)
idx3 = KLIEPInference.trimap(4, 5)
idx4 = KLIEPInference.trimap(4, 6)

# m = 25
γx = chain(25, 25, 0.1, 0.3, 0)
γy = deepcopy(γx)
γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
γy[idx1] *= -1.0
γy[idx2] *= -1.0
γy[idx3] -= 0.1
γy[idx4] -= 0.1
@save "graphs/params_exp1_chain2_25.jld" γx γy

# m = 50
γx = chain(50, 50, 0.1, 0.3, 0)
γy = deepcopy(γx)
γy[idx0] = abs(γy[idx0]) > 0.2 ? γy[idx0] - sign(γy[idx0])*0.1 : γy[idx0] + sign(γy[idx0])*0.1
γy[idx1] *= -1.0
γy[idx2] *= -1.0
γy[idx3] -= 0.1
γy[idx4] -= 0.1
@save "graphs/params_exp1_chain2_50.jld" γx γy
