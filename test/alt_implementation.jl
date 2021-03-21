module Alt

using KLIEPInference
using LinearAlgebra
import JuMP, SCS, MathOptInterface, MosekTools

export
  SCS_KLIEP,
  Mosek_KLIEP

struct SCS_KLIEP <: KLIEPSolver end
struct Mosek_KLIEP <: KLIEPSolver end

KLIEP(Ψx, Ψy, ::SCS_KLIEP) = 
  _KLIEP(Ψx, Ψy, JuMP.optimizer_with_attributes(SCS.Optimizer, "verbose" => 1))
KLIEP(Ψx, Ψy, ::Mosek_KLIEP) = 
  _KLIEP(Ψx, Ψy, JuMP.optimizer_with_attributes(MosekTools.Mosek.Optimizer, "QUIET" => false))

spKLIEP(Ψx, Ψy, λ, ::SCS_KLIEP) =
  _spKLIEP!(Vector{Float64}(undef, size(Ψx, 1)), Ψx, Ψy, λ, JuMP.optimizer_with_attributes(SCS.Optimizer, "verbose" => 1))
spKLIEP(Ψx, Ψy, λ, ::Mosek_KLIEP) =
  _spKLIEP!(Vector{Float64}(undef, size(Ψx, 1)), Ψx, Ψy, λ, JuMP.optimizer_with_attributes(MosekTools.Mosek.Optimizer, "QUIET" => false))


### solvers

function _KLIEP(Ψx, Ψy, solver)
  m, nx = size(Ψx)
  ny = size(Ψy, 2)
  @assert m == size(Ψy, 1)

  lny = log(ny)

  problem = JuMP.Model(solver)
  JuMP.@variable(problem, Θ[1:m])
  JuMP.@variable(problem, t)
  JuMP.@variable(problem, u[1:ny])

  JuMP.@constraint(problem, sum(u) <= 1.)
  for j=1:ny
      JuMP.@constraint(problem, [dot(Θ, Ψy[:, j]) - t - lny, 1., u[j]] in MathOptInterface.ExponentialCone())
  end

  JuMP.@objective(problem, Min, -sum(Ψx'*Θ) / nx + t)
  JuMP.optimize!(problem)

  JuMP.value.(Θ)
end

function _spKLIEP!(θhat, Ψx, Ψy, λ, solver)

  m, nx = size(Ψx)
  ny = size(Ψy, 2)
  @assert m == size(Ψy, 1) == length(θhat)

  lny = log(ny)

  problem = JuMP.Model(solver)
  JuMP.@variable(problem, x[1:m])
  JuMP.@variable(problem, t)
  JuMP.@variable(problem, u[1:ny])
  JuMP.@variable(problem, b[1:m])

  # log-sum-exp constraints
  JuMP.@constraint(problem, sum(u) <= 1.)
  for j=1:ny
      JuMP.@constraint(problem, [dot(x, Ψy[:, j]) - t - lny, 1., u[j]] in MathOptInterface.ExponentialCone())
  end

  # l1 constraints
  JuMP.@constraint(problem,   x .<= b )
  JuMP.@constraint(problem,  -x .<= b )

  JuMP.@objective(problem, Min, -sum(Ψx'*x) / nx + t + λ*sum(b))
  JuMP.optimize!(problem)

  θhat .= JuMP.value.(x)
end

end
