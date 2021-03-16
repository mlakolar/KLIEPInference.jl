m = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
numChanges = parse(Int,ARGS[3])
lbInd = parse(Int,ARGS[4])
rep = parse(Int, ARGS[5])

if isfile("/scratch/midway2/byolkim/exper4/res_$(m)_$(sgn)_$(numChanges)_$(lbInd)_$(rep).jld")
    println("the file already exists!")
    exit()
end

using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, StatsBase, JLD

println("importing parameters from params_exp4_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld...")
file = jldopen("params_exp4_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld", "r")
γx = read(file, "γx")
γy = read(file, "γy")
close(file)

p = length(γx)
nx = 500
ny = 500

println("generating samples...")
Random.seed!(123 + rep)
spl = IsingSampler(γx; thin=2000)
X = rand(spl, nx)
spl = IsingSampler(γy; thin=2000)
Y = rand(spl, ny)

Ψx = Ψising(X)
Ψy = Ψising(Y)

println("step 1")
λ1 = 1.01 * quantile(Normal(), 1. - 0.05 / p)
θ = spKLIEP(Ψx, Ψy, λ1, CD_KLIEP(); loadings=true)

println("step 2")
λ2 = sqrt(2. * log(p) / ny)
H = KLIEP_Hessian(spzeros(Float64, p), Ψy)
Hinv = Vector{SparseIterate{Float64}}(undef, p)
for k = 1:p
    ω = Hinv_row(H, k, λ2)

    supp = KLIEPInference._find_supp(k, ω)
    h = view(H, supp, supp)
    δ = (supp .== k)
    ω[supp] = h\δ

    Hinv[k] = ω
end

println("step 3 + bootstrap...")
boot1, boot2 = boot_SparKLIE(Ψx, Ψy, θ, Hinv)

println("step 3 + bootstrap completed, testing at α = 0.05...")
CI = simulCI(boot1, 0.05)
T1 = all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

CI = simulCI(boot2, 0.05)
T2 = all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

CI = simulCIstudentized(boot1, 0.05)
W1 = all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

CI = simulCIstudentized(boot2, 0.05)
W2 = all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

@save "/scratch/midway2/byolkim/exper4/res_$(m)_$(sgn)_$(numChanges)_$(lbInd)_$(rep).jld" T1 T2 W1 W2
