using KLIEPInference
using Distributions
using LinearAlgebra
using SparseArrays
using Random
using JLD

fname = ARGS[1]
nx = parse(Int,ARGS[2])
ny = parse(Int,ARGS[3])
rep = parse(Int, ARGS[4])

file = jldopen(fname, "r")
p = read(file, "p")
sgn = read(file, "sgn")
θx = read(file, "θx")
θy = read(file, "θy")
close(file)

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
# step 1
#
###########################
@show "step 1"
m = div(p * (p - 1),  2)
@show λ1 = 1.2 * sqrt(2. * log(m) / nx)
θhat = spKLIEP(Ψx, Ψy, λ1, CD_KLIEP())
spKLIEP_refit!(θhat, Ψx, Ψy)
@show "step 1 done"

###########################
#
# step 2
#
###########################
@show "step 2"
@show λ2 = 1.2 * sqrt(2. * log(p) / ny)
H = KLIEP_Hessian(θhat, Ψy)
Hinv = Vector{SparseVector{Float64,Int64}}(undef, m)
for row=1:m
    ω = Hinv_row(H, row, λ2)
    Hinv[row] = convert(SparseVector, ω)
end
@show "step 2 done"


###########################
#
# step 3
#
###########################

res = boot_spKLIEP(Ψx, Ψy, θhat, Hinv; bootSamples=300)

@save "/scratch/midway2/mkolar/KLIEP/exp1/res_$(p)_$(sgn)_$(nx)_$(ny)_$(rep).jld" θhat Hinv res
