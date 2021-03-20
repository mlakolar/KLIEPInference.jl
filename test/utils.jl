module UtilTest

using Test
using Statistics
import KLIEPInference: trimap, itrimap, pack, unpack, Ψising

@testset "trimap" begin
    @test trimap(2, 1) == 1
	@test trimap(1, 2) == 1
	@test trimap(4, 3) == 6

	@test itrimap(1) == CartesianIndex(2, 1)
	@test itrimap(6) == CartesianIndex(4, 3)
end

@testset "pack" begin
	A = [1. 2.; 2. 3.]
	θ = [2.]
	@test pack(A) == θ
	@test unpack(θ) == [0. 2.; 2. 0]

	A = [1. 2. 3.; 2. 2. 0.5; 3. 0.5 3.]
	θ = [2., 3., 0.5]
	@test pack(A) == θ
	@test unpack(θ) == [0. 2. 3.; 2. 0. 0.5; 3. 0.5 0.]
end

@testset "Ising" begin
	X = Array{Bool}([1 1; 0 1; 0 0; 1 0])
	Ψ = [-1. 1; -1 -1; 1 -1; 1 -1; -1 -1; -1 1]	
	@test (@test_logs (:warn, "a row of Ψ has zero variance") Ψising(X)) == Ψ


end

end



