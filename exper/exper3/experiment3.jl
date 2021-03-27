m = parse(Int, ARGS[1])
sgn = parse(Int, ARGS[2])
rep = parse(Int, ARGS[3])

scratch_dir = ARGS[4]   # e.g. "/scratch/midway2/byolkim/exper3"

if isfile("$(scratch_dir)/res_$(m)_$(sgn)_$(rep).jld")
    println("the file already exists!")
    exit()
end

using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, StatsBase, JLD

println("importing parameters from params_exp3_$(m)_$(sgn).jld")
file = jldopen("params_exp3_$(m)_$(sgn).jld", "r")
γx = read(file, "γx")
γy = read(file, "γy")
close(file)

p = length(γx)
nx = 500
ny = 500

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
Hinv = Vector{SparseIterate{Float64}}(undef, p)
for k = 1:p
    ω = Hinv_row(H, k, λ2)

    supp = KLIEPInference._find_supp(k, ω)
    h = view(H, supp, supp)
    δ = (supp .== k)
    ω[supp] = h\δ

    Hinv[k] = ω
end

println("step 3 + bootstrap")
boot1, boot2 = boot_SparKLIE(Ψx, Ψy, θ, Hinv; bootSamples=1000)

println("saving results to $(scratch_dir)/res_$(m)_$(sgn)_$(rep).jld")
@save "$(scratch_dir)/res_$(m)_$(sgn)_$(rep).jld" boot1 boot2

println("done!")
