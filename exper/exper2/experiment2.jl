set = ARGS[1]
δ = parse(Float64, ARGS[2])
batch = parse(Int, ARGS[3])
#rep = parse(Int, ARGS[3])
scratch_dir = ARGS[4]   # e.g. "/scratch/midway2/byolkim/exper2"

if isfile("$(scratch_dir)/res_$(set)_$(δ)_$(batch).jld")
    println("the file already exists!")
    exit()
end

using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, StatsBase, JLD

function generate_data(γx, γy, rep)
    nx = 150
    ny = 300

    println("Repetition $(rep)")
    println("=================")
    println("generating samples")
    Random.seed!(123 + rep)
    spl = IsingSampler(γx; thin=2000)
    X = rand(spl, nx)
    spl = IsingSampler(γy; thin=2000)
    Y = rand(spl, ny)

    Ψx = Ψising(X)
    Ψy = Ψising(Y)

    Ψx, Ψy
end

function test_exper2(Ψx, Ψy, p)
    nx = 150
    ny = 300
    idx = KLIEPInference.trimap(5, 6)    

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

    println("done!")
    T1, T2
end


function exper2(set, δ, batch, scratch_dir)    
    println("importing parameters from params_exp2_$(set)_$(δ).jld")
    file = jldopen("graphs/params_exp2_$(set)_$(δ).jld", "r")
    γx = read(file, "γx")
    γy = read(file, "γy")
    close(file)
    p = length(γx)

    T1 = zeros(100)
    T2 = zeros(100)

    start_rep = (batch-1) * 100
    for rep=1:100
        Ψx, Ψy = generate_data(γx, γy, start_rep + rep)
        T1[rep], T2[rep] = test_exper2(Ψx, Ψy, p)
    end

    println("saving results to $(scratch_dir)/res_$(set)_$(δ)_$(batch).jld")
    @save "$(scratch_dir)/res_$(set)_$(δ)_$(batch).jld" T1 T2
    nothing    
end

exper2(set, δ, batch, scratch_dir)





