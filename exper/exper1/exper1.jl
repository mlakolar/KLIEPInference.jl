module Exper1

using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics, Random
using Distributions, JLD
using KLIEPInference: trimap


function experiment1(Θx, Θy, idx, nx, ny, nrep)
    p = length(Θx)

    # for both steps, we use autoscaling procedures with canonical choices of λ
    λ1 = 1.01 * quantile(Normal(), 1. - 0.05 / p)
    λ2 = sqrt(2. * log(p) / ny)

    res = zeros(Float64, nrep, 4, 2)
    for rep = 1:nrep
        if mod(rep, div(nrep, 10)) === 0
            println("$(rep) / $(nrep)")
        end

        # generate samples
        spl = IsingSampler(Θx; thin=2000)
        X = rand(spl, nx)
        spl = IsingSampler(Θy; thin=2000)
        Y = rand(spl, ny)

        Ψx = Ψising(X)
        Ψy = Ψising(Y)

        # oracle estimate
        supp = findall(!iszero, Θx - Θy)
        θ = KLIEP(Ψx[supp, :], Ψy[supp, :], CD_KLIEP())

        ω = KLIEP_Hessian(θ, Ψy[supp, :])\(supp .=== idx)
        σ = stderr_SparKLIE(Ψx[supp, :], Ψy[supp, :], θ, ω)

        res[rep, 1, 1] = θ[findfirst(isequal(idx), supp)] - Θx[idx] + Θy[idx]
        res[rep, 1, 2] = res[rep, 1, 1] / σ

        # step 1
        θ = spKLIEP(Ψx, Ψy, λ1, CD_KLIEP(); loadings=true)

        # naïve re-fitted estimate
        supp = KLIEPInference._find_supp(idx, θ)
        θ[supp] = KLIEP(Ψx[supp, :], Ψy[supp, :], CD_KLIEP())

        ω = KLIEP_Hessian(θ[supp], Ψy[supp, :])\(supp .=== idx)
        σ = stderr_SparKLIE(Ψx[supp, :], Ψy[supp, :], θ[supp], ω)

        res[rep, 2, 1] = θ[idx] - Θx[idx] + Θy[idx]
        res[rep, 2, 2] = res[rep, 2, 1] / σ

        # step 2
        H = KLIEP_Hessian(θ, Ψy)
        ω = Hinv_row(H, idx, λ2)
        supp = KLIEPInference._find_supp(idx, ω)
        ω[supp] = view(H, supp, supp)\(supp .=== idx)

        # SparKLIE+1
        θ1 = KLIEPInference._debias1(Ψx, Ψy, θ, ω, idx; refit=false)
        σ = stderr_SparKLIE(Ψx, Ψy, θ, ω)

        res[rep, 3, 1] = θ1 - Θx[idx] + Θy[idx]
        res[rep, 3, 2] = res[rep, 3, 1] / σ

        # SparKLIE+2
        θ2 = KLIEPInference._debias2(Ψx, Ψy, θ, ω, idx)

        supp = KLIEPInference._find_supp(idx, θ, ω)
        ω = KLIEP_Hessian(θ[supp], Ψy[supp, :])\(supp .=== idx)
        σ = stderr_SparKLIE(Ψx[supp, :], Ψy[supp, :], θ[supp], ω)

        res[rep, 4, 1] = θ2 - Θx[idx] + Θy[idx]
        res[rep, 4, 2] = res[rep, 4, 1] / σ
    end
    res
end

function print_res(res)
    println("avg bias: $(round.(mean(res[:, :, 1], dims=1), digits=5))")
    println("coverage: $(mean(abs.(res[:, :, 2]) .< quantile(Normal(), 0.975), dims=1))")
end

function chain1_25()
    Random.seed!(325)
    m = 25 
    γx = chain(m; lenC=m, lb=0.1, ub=0.2, sgn=0)    
    Δ = zeros(length(γx)) 
    idx0 = trimap(5, 6)
    idx1 = trimap(4, 5)
    idx2 = trimap(6, 7)
    idx3 = trimap(4, 6)
    idx4 = trimap(5, 7)   

    Δ[idx0] = -0.2
    Δ[idx1] =  0.4
    Δ[idx2] = -0.4
    Δ[idx3] = 0.2
    Δ[idx4] = 0.2

    γy = γx - Δ
    return γx, γy
end

function run_chain1_25()
    γx, γy = chain1_25()
    nx = 150
    ny = 300 
    res = experiment1(γx, γy, trimap(5, 6), nx, ny, 1000)
    @save "new_res_chain1_25.jld" res 
    print_res(res)
end

run_chain1_25()

end