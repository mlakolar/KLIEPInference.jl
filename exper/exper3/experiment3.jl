using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, StatsBase, JLD

m = parse(Int, ARGS[1])
sgn = parse(Int, ARGS[2])
rep = parse(Int, ARGS[3])

println("importing parameters from params_exp3_$(m)_$(sgn).jld...")
file = jldopen("params_exp3_$(m)_$(sgn).jld", "r")
θx = read(file, "θx")
θy = read(file, "θy")
close(file)

p = length(θx)
nx = 500
ny = 500

println("generating samples...")
Random.seed!(123 + rep)
spl = IsingSampler(θx; thin=2000)
X = rand(spl, nx)
spl = IsingSampler(θy; thin=2000)
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

println("step 3 + bootstrap completed, computing coverage...")
α = 0.05:0.05:0.95
nα = length(α)
res = zeros(Bool, 2, nα)
for i in 1:nα
    CI = simulCI(boot1, α[i])
    res[1, i] = all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

    CI = simulCI(boot2, α[i])
    res[2, i] = all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0
end

println("saving results to /scratch/midway2/byolkim/exper3/res_$(m)_$(sgn)_$(rep).jld")
@save "/scratch/midway2/byolkim/exper3/res_$(m)_$(sgn)_$(rep).jld" res

println("done!")
