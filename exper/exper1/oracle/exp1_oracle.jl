using KLIEPInference
using Distributions
using LinearAlgebra
using SparseArrays
using Random
using JLD

p = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
nx = parse(Int,ARGS[3])
ny = parse(Int,ARGS[4])
rep = parse(Int, ARGS[5])

file = jldopen("params_exp1_$(p)_$(sgn).jld", "r")
θx = read(file, "θx")
θy = read(file, "θy")
close(file)

Δ = θx - θy

# generate data

Random.seed!(123 + rep)
spl = IsingSampler(θx; thin=2000)
X = rand(spl, nx)
spl = IsingSampler(θy; thin=2000)
Y = rand(spl, ny)
Ψx = Ψising(X)
Ψy = Ψising(Y)


###########################
#
# bootstrap step
#
###########################

S = findall(!iszero, Δ)
res = boot_oracleKLIEP(Ψx, Ψy, S; bootSamples=300)

@save "/scratch/midway2/mkolar/KLIEP/exp1/oracle/res_$(p)_$(sgn)_$(nx)_$(ny)_$(rep).jld" res
