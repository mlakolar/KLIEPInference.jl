module SolverTest

using Test, Random, Statistics, SparseArrays
import KLIEPInference
import KLIEPInference: chain, IsingSampler, Ψising, KLIEP, CD_KLIEP          
import KLIEPInference: spKLIEP       
using Distributions

include("alt_implementation.jl")

import .Alt

@testset "KLIEP low dim" begin
    Random.seed!(3)
    θx = chain(6, 3, 0.2, 0.4, 0)
    θy = copy(θx)

    nx = 500
    ny = 500    
    X = rand(IsingSampler(θx; thin=1000), nx)    
    Y = rand(IsingSampler(θy; thin=1000), ny)
    Ψx = Ψising(X)
    Ψy = Ψising(Y)

    scs_out = Alt.KLIEP(Ψx, Ψy, Alt.SCS_KLIEP())
    cd_out = KLIEP(Ψx, Ψy, CD_KLIEP())
    @test scs_out ≈ cd_out atol=0.001

    # mosek_out = Alt.KLIEP(Ψx, Ψy, Alt.Mosek_KLIEP())
    # @test scs_out ≈ mosek_out atol=0.01
end


@testset "spKLIEP delta == 0" begin
    Random.seed!(4)
    θx = chain(30, 10, 0.1, 0.3, 0)
    θy = copy(θx)

    nx = 300
    ny = 300    
    X = rand(IsingSampler(θx; thin=1000), nx)    
    Y = rand(IsingSampler(θy; thin=1000), ny)
    Ψx = Ψising(X)
    Ψy = Ψising(Y)

    p = length(θx)
    λ = 1.01 * quantile(Normal(), 1. - 0.05 / p) / sqrt(nx)

    θ = spKLIEP(Ψx, Ψy, λ, CD_KLIEP(); loadings=false)

    scs_out = Alt.spKLIEP(Ψx, Ψy, λ, Alt.SCS_KLIEP())    
    @test scs_out ≈ θ atol=0.001    

    # mosek_out = Alt.spKLIEP(Ψx, Ψy, λ, Alt.Mosek_KLIEP())
    # @test mosek_out ≈ θ atol=0.001    
end


@testset "spKLIEP delta != 0" begin
    Random.seed!(5)
    θx = chain(30, 10, 0.1, 0.3, 0)
    θy = copy(θx)
    p = length(θx)  
    Δ = sprand(p, 0.1)

    du = Uniform(0.1, 0.2)
    I, V = findnz(Δ)
    for i in I        
        nv = rand(du)
        nv *= 2. * rand(Bernoulli()) - 1.
        θy[i] += nv
    end
    Δ = sparse(θx - θy)

    nx = 500
    ny = 500    
    X = rand(IsingSampler(θx; thin=1000), nx)    
    Y = rand(IsingSampler(θy; thin=1000), ny)
    Ψx = Ψising(X)
    Ψy = Ψising(Y)
    
    λ = 1.01 * quantile(Normal(), 1. - 0.05 / p) / sqrt(nx)

    θ = spKLIEP(Ψx, Ψy, λ, CD_KLIEP(); loadings=false)

    scs_out = Alt.spKLIEP(Ψx, Ψy, λ, Alt.SCS_KLIEP())    
    @test scs_out ≈ θ atol=0.001    

    # mosek_out = Alt.spKLIEP(Ψx, Ψy, λ, Alt.Mosek_KLIEP())
    # @test mosek_out ≈ θ atol=0.001    
end



end