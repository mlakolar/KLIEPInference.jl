module SolverTest

using Test, Random
import KLIEPInference
import KLIEPInference: chain, IsingSampler, Ψising, KLIEP, CD_KLIEP          

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
    # mosek_out = Alt.KLIEP(Ψx, Ψy, Alt.Mosek_KLIEP())
    cd_out = KLIEP(Ψx, Ψy, CD_KLIEP())

    # @test scs_out ≈ mosek_out atol=0.01
    @test scs_out ≈ cd_out atol=0.01
end


end