set = "set1"
δ = 0.0
scratch_dir = "tmp"

if isfile("$(scratch_dir)/res_$(set)_$(δ)_$(rep).jld")
    println("the file already exists!")
    exit()
end

using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, StatsBase, JLD

println("importing parameters from params_exp2_$(set)_$(δ).jld")
file = jldopen("graphs/params_exp2_$(set)_$(δ).jld", "r")
γx = read(file, "γx")
γy = read(file, "γy")
close(file)

idx = KLIEPInference.trimap(5, 6)

p = length(γx)
nx = 150
ny = 300

println("generating samples")
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
H = KLIEP_Hessian(θ, Ψy)
ω = Hinv_row(H, idx, λ2)
supp = KLIEPInference._find_supp(idx, ω)
ω[supp] = view(H, supp, supp)\(supp .=== idx)

println("step 3")
# SparKLIE+1
θ1 = KLIEPInference._debias1(Ψx, Ψy, θ, ω, idx)
σ = stderr_SparKLIE(Ψx, Ψy, θ, ω)

T1 = abs(θ1 / σ) > quantile(Normal(), 0.975) ? 1 : 0

# SparKLIE+2
θ2 = KLIEPInference._debias2(Ψx, Ψy, θ, ω, idx)
supp = KLIEPInference._find_supp(idx, θ, ω)
ω = KLIEP_Hessian(θ[supp], Ψy[supp, :])\(supp .=== idx)
σ = stderr_SparKLIE(Ψx[supp, :], Ψy[supp, :], θ[supp], ω)

T2 = abs(θ2 / σ) > quantile(Normal(), 0.975) ? 1 : 0

println("saving results to $(scratch_dir)/res_$(set)_$(δ)_$(rep).jld")
@save "$(scratch_dir)/res_$(set)_$(δ)_$(rep).jld" T1 T2

println("done!")
