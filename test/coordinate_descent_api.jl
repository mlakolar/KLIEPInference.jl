module CoordinateDescentAPITest
using Test, CoordinateDescent, ProximalBase, KLIEPInference, LinearAlgebra, Statistics

# Independent binary features give a diagonal Hessian at zero. Their exact
# tilting parameter is atanh(mu), so no reference optimizer is needed.
@testset "CoordinateDescent 0.8 interface" begin
  Ψy = [-1. -1. 1. 1.; -1. 1. -1. 1.]
  Ψx = [0.2 0.2 0.2 0.2; -0.3 -0.3 -0.3 -0.3]
  x = SparseIterate(2)
  @test KLIEP!(x, Ψx, Ψy) === x
  @test Vector(x) ≈ atanh.([0.2, -0.3]) atol=1e-6

  x = SparseIterate(2)
  @test spKLIEP!(x, Ψx, Ψy, 0.1) === x
  @test Vector(x) ≈ atanh.([0.1, -0.2]) atol=1e-6
  f = KLIEPInference.CDKLIEPLoss(Ψx, Ψy)
  @test CoordinateDescent.num_coordinates(f) == 2
  @test spKLIEP!(x, f, NormL1(0.1)) === x
  @test Vector(x) ≈ atanh.([0.1, -0.2]) atol=1e-6
  @test Hinv_row(Matrix(1.0I, 3, 3), 2, 0.1) ≈ [0., 1., 0.]
end
end
