module UtilTest

using Test
using Statistics, LinearAlgebra
import KLIEPInference: trimap, itrimap, pack, unpack, Ψising, rhat, KLIEP_Hessian, _find_supp
import KLIEPInference: chain, IsingSampler
import Random: rand

@testset "trimap" begin
    @test trimap(1, 2) == 1
	@test trimap(3, 4) == 6

	@test itrimap(1) == CartesianIndex(1, 2)
	@test itrimap(6) == CartesianIndex(3, 4)
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
	@test (@test_logs (:warn, "a row of Ψ has zero variance") match_mode=:any Ψising(X)) == Ψ

	for rep=1:50
		n = 200
		p = 3
		m = 3
		θ = chain(m; lenC=3, lb=0.1, ub=0.2, sgn=0)
		spl = IsingSampler(θ; burn=5000, thin=5000)
		X = rand(spl, n)
		Ψ = Ψising(X)

		# check rhat
		_rhat = zeros(n)
		_mr = 0.
		for i = 1:n
			_rhat[i] = exp( dot(Ψ[:, i], θ) )
			_mr += _rhat[i]
		end
		_mr /= n
		for i = 1:n
			_rhat[i] /= _mr
		end
		@test rhat(θ, Ψ) ≈ _rhat atol=1e-5

		# check Hessian
		_H = zeros(3, 3)
		_hμ = zeros(3)
		for j=1:p
			for i=1:n		
				_hμ[j] += Ψ[j, i] * _rhat[i]
			end
			_hμ[j] /= n
		end
		for j=1:p
			for k=1:p
				for i=1:n
					_H[j, k] += Ψ[j, i] * Ψ[k, i] * _rhat[i]
				end
				_H[j, k] /= n
			end
		end
		_H -= _hμ * transpose(_hμ)
		@test KLIEP_Hessian(θ, Ψ) ≈ _H atol=1e-5 
	end
end

@testset "_find_supp" begin
	x = [0., 1., 0., 3., 0.]
	@test _find_supp(1, x) == [2, 4, 1]
	@test _find_supp(2, x) == [4, 2]	
	@test _find_supp(3, x) == [2, 4, 3]
	@test _find_supp(4, x) == [2, 4]
	@test _find_supp(5, x) == [2, 4, 5]

	x = [0., 0., 0., 3., 0.]
	y = [0., 1., 1., 3., 0.]
	@test _find_supp(1, x, y) == [4, 2, 3, 1]
	@test _find_supp(2, x, y) == [4, 3, 2]
	@test _find_supp(3, x, y) == [4, 2, 3]
	@test _find_supp(4, x, y) == [3, 2, 4]
	@test _find_supp(5, x, y) == [4, 2, 3, 5]
end


end



