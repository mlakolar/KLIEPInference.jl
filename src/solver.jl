abstract type KLIEPSolver end

struct CD_KLIEP <: KLIEPSolver end
struct SCS_KLIEP <: KLIEPSolver end
struct Mosek_KLIEP <: KLIEPSolver end

### interface

KLIEP(Ψx, Ψy, ::SCS_KLIEP) = _KLIEP(Ψx, Ψy, JuMP.with_optimizer(SCS.Optimizer, verbose=1))
KLIEP(Ψx, Ψy, ::Mosek_KLIEP) = _KLIEP(Ψx, Ψy, JuMP.with_optimizer(MathOptInterfaceMosek.MosekOptimizer, QUIET=false))

KLIEP(Ψx, Ψy, ::CD_KLIEP) = KLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy)

function KLIEP!(x::SparseIterate, Ψx, Ψy)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxZero()

    coordinateDescent!(x, f, g)
end





spKLIEP(Ψx, Ψy, λ, ::SCS_KLIEP) =
  _spKLIEP!(Vector{Float64}(undef, size(Ψx, 1)), Ψx, Ψy, λ, JuMP.with_optimizer(SCS.Optimizer, verbose=1))

function spKLIEP(Ψx, Ψy, λ, ::CD_KLIEP)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxL1(λ)

    convert(Vector, coordinateDescent!(SparseIterate(f.p), f, g))
end





function Hinv_row(H, row, λ0)
    m = size(H, 1)
    e = zeros(m)
    e[row] = -1.0

    x = SparseIterate(m)
    x[row] = 1.
    f = CDQuadraticLoss(H, e)

    # init sigma
    σ = sqrt( dot(x, H * x) )

    for iter=1:10
        coordinateDescent!(x, f, ProxL1(λ0 * σ))
        σnew = sqrt( dot(x, H * x) )

        if abs(σnew - σ) / σ < 1e-2
          break
        end
        σ = σnew
    end

    convert(Vector, x)
end


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
        JuMP.@constraint(problem, [dot(Θ, Ψy[:, j]) - t - lny, 1., u[j]] in MOI.ExponentialCone())
    end

    JuMP.@objective(problem, Min, -sum(Ψx'*Θ) / nx + t)
    JuMP.optimize!(problem)

    JuMP.result_value.(Θ)
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
        JuMP.@constraint(problem, [dot(x, Ψy[:, j]) - t - lny, 1., u[j]] in MOI.ExponentialCone())
    end

    # l1 constraints
    JuMP.@constraint(problem,   x .<= b )
    JuMP.@constraint(problem,  -x .<= b )

    JuMP.@objective(problem, Min, -sum(Ψx'*x) / nx + t + λ*sum(b))
    JuMP.optimize!(problem)

    θhat .= JuMP.result_value.(x)

end







####################################
#
# loss En[-θ'Ψx(i)] + log( En[exp(θ'Ψy(i))] )
#
####################################
struct CDKLIEPLoss <: CoordinateDifferentiableFunction
  μx::Vector{Float64}   # mean(Ψx, dims = 2)
  Ψy::Matrix{Float64}
  r::Vector{Float64}    # ny dim vector --- stores θ'Ψy(i)
  p::Int64
end

function CDKLIEPLoss(Ψx::Matrix{Float64}, Ψy::Matrix{Float64})
  (p = size(Ψx, 1)) == size(Ψy, 1) || throw(DimensionMismatch())

  CDKLIEPLoss(vec(mean(Ψx, dims=2)), Ψy, zeros(size(Ψy, 2)), p)
end

CoordinateDescent.numCoordinates(f::CDKLIEPLoss) = f.p

function CoordinateDescent.initialize!(f::CDKLIEPLoss, x::SparseIterate)

  mul!(f.r, transpose(f.Ψy), x)

  nothing
end

function CoordinateDescent.gradient(
  f::CDKLIEPLoss,
  x::SparseIterate,
  j::Int64)

  μx, Ψy, r = f.μx, f.Ψy, f.r

  -μx[j] + mean( exp.(r) .* Ψy[j, :] ) / mean(exp, r)
end


function CoordinateDescent.descendCoordinate!(
  f::CDKLIEPLoss,
  g::Union{ProxL1, ProxZero},
  x::SparseIterate,
  j::Int64)

  μx, Ψy, r = f.μx, f.Ψy, f.r

  mean_exp_r = mean(exp, r)
  mean_exp_r_Ψy = mean( exp.(r) .* Ψy[j, :] )
  mean_exp_r_Ψy2 = mean( exp.(r) .* Ψy[j, :].^2. )

  @inbounds grad = -μx[j] + mean_exp_r_Ψy / mean_exp_r
  H = mean_exp_r_Ψy2 / mean_exp_r - mean_exp_r_Ψy^2 / mean_exp_r^2.

  @inbounds oldVal = x[j]
  @inbounds x[j] -= grad / H
  newVal = cdprox!(g, x, j, 1. / H)
  h = newVal - oldVal

  # update internals
  for i=1:length(r)
      @inbounds r[i] += h*Ψy[j, i]
  end
  h
end
