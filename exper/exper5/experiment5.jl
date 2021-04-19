graph = ARGS[1]
idx = parse(Int, ARGS[2])
nx = parse(Int, ARGS[3])
ny = parse(Int, ARGS[4])
rep = parse(Int, ARGS[5])

scratch_dir = ARGS[6]

if isfile("$(scratch_dir)/res_$(graph)_$(rep).jld")
    println("the file already exists!")
    exit()
end

using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, StatsBase, JLD

include("SymKLIEP.jl")

println("importing parameters from $(graph).jld")
file = jldopen("graphs/$(graph).jld", "r")
γx = pack(read(file, "Θx"))
γy = pack(read(file, "Θy"))
close(file)

p = length(γx)

println("generating samples")
Random.seed!(123 + rep)
spl = IsingSampler(γx; thin=2000)
X = rand(spl, nx)
spl = IsingSampler(γy; thin=2000)
Y = rand(spl, ny)

Ψx = Ψising(X)
Ψy = Ψising(Y)

res = zeros(3, 2, 3, 5)

println("original procedure with ny > nx ...")
@show λ1 = exp.(range(log(sqrt(16. * log(p) / nx)), stop=log(sqrt(2. * log(p) / nx)), length=5))
λ2 = sqrt(2. * log(p) / ny)
for t = 1:5
    # step 1
    θ = spKLIEP(Ψx, Ψy, λ1[t], CD_KLIEP(); loadings=false)

    # step 2
    H = KLIEP_Hessian(θ, Ψy)
    ω = Hinv_row(H, idx, λ2)
    supp = KLIEPInference._find_supp(idx, ω)
    ω[supp] = view(H, supp, supp)\(supp .=== idx)

    # SparKLIE+1
    θ1 = KLIEPInference._debias1(Ψx, Ψy, θ, ω, idx)
    σ = stderr_SparKLIE(Ψx, Ψy, θ, ω)

    res[1, 1, 1, t] = θ1 - γx[idx] + γy[idx]
    res[1, 1, 2, t] = σ
    res[1, 1, 3, t] = (θ1 - γx[idx] + γy[idx]) / σ

    # SparKLIE+2
    θ2 = KLIEPInference._debias2(Ψx, Ψy, θ, ω, idx)

    supp = KLIEPInference._find_supp(idx, θ, ω)
    ω = KLIEP_Hessian(θ[supp], Ψy[supp, :])\(supp .=== idx)
    σ = stderr_SparKLIE(Ψx[supp, :], Ψy[supp, :], θ[supp], ω)

    res[1, 2, 1, t] = θ2 - γx[idx] + γy[idx]
    res[1, 2, 2, t] = σ
    res[1, 2, 3, t] = (θ2 - γx[idx] + γy[idx]) / σ
end
println("... done!")

println("reverse KL ...")
@show λ1 = exp.(range(log(sqrt(16. * log(p) / nx)), stop=log(sqrt(2. * log(p) / nx)), length=5))
λ2 = sqrt(2. * log(p) / nx)
for t = 1:5
    # step 1
    θ = spKLIEP(Ψy, Ψx, λ1[t], CD_KLIEP(); loadings=false)

    # step 2
    H = KLIEP_Hessian(θ, Ψx)
    ω = Hinv_row(H, idx, λ2)
    supp = KLIEPInference._find_supp(idx, ω)
    ω[supp] = view(H, supp, supp)\(supp .=== idx)

    # SparKLIE+1
    θ1 = KLIEPInference._debias1(Ψy, Ψx, θ, ω, idx)
    σ = stderr_SparKLIE(Ψy, Ψx, θ, ω)

    res[2, 1, 1, t] = -θ1 + γx[idx] - γy[idx]
    res[2, 1, 2, t] = σ
    res[2, 1, 3, t] = (-θ1 + γx[idx] - γy[idx]) / σ

    # SparKLIE+2
    θ2 = KLIEPInference._debias2(Ψy, Ψx, θ, ω, idx)

    supp = KLIEPInference._find_supp(idx, θ, ω)
    ω = KLIEP_Hessian(θ[supp], Ψx[supp, :])\(supp .=== idx)
    σ = stderr_SparKLIE(Ψy[supp, :], Ψx[supp, :], θ[supp], ω)

    res[2, 2, 1, t] = -θ2 + γx[idx] - γy[idx]
    res[2, 2, 2, t] = σ
    res[2, 2, 3, t] = (-θ2 + γx[idx] - γy[idx]) / σ
end
println("... done!")

println("symmetric KL ...")
@show λ1 = exp.(range(log(sqrt(16. * log(p) / nx)), stop=log(sqrt(2. * log(p) / nx)), length=5))
λ2 = .5 * (sqrt(2. * log(p) / nx) + sqrt(2. * log(p) / ny))
for t = 1:5
    # step 1
    θ = spSymKLIEP(Ψx, Ψy, λ1[t], CD_SymKLIEP())

    # step 2
    H = SymKLIEP_Hessian(θ, Ψx, Ψy)
    ω = Hinv_row(H, idx, λ2)
    supp = KLIEPInference._find_supp(idx, ω)
    ω[supp] = view(H, supp, supp)\(supp .=== idx)

    # SparKLIE+1
    θ1 = SymKLIEP_debias1(Ψx, Ψy, θ, ω, idx)
    σ = SymKLIEP_stderr(Ψx, Ψy, θ, ω)

    res[3, 1, 1, t] = θ1 - γx[idx] + γy[idx]
    res[3, 1, 2, t] = σ
    res[3, 1, 3, t] = (θ1 - γx[idx] + γy[idx]) / σ

    # SparKLIE+2
    θ2 = SymKLIEP_debias2(Ψx, Ψy, θ, ω, idx)

    supp = KLIEPInference._find_supp(idx, θ, ω)
    ω = SymKLIEP_Hessian(θ[supp], Ψx[supp, :], Ψy[supp, :])\(supp .=== idx)
    σ = SymKLIEP_stderr(Ψx[supp, :], Ψy[supp, :], θ[supp], ω)

    res[3, 2, 1, t] = θ2 - γx[idx] + γy[idx]
    res[3, 2, 2, t] = σ
    res[3, 2, 3, t] = (θ2 - γx[idx] + γy[idx]) / σ
end
println("... done!")

println("saving results to $(scratch_dir)/res_$(graph)_$(rep).jld ...")
@save "$(scratch_dir)/res_$(graph)_$(rep).jld" res
println("... done!")
