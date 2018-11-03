
function Ψising(X)
    m, n = size(X)

    p = div(m * (m - 1), 2)
    out = zeros(Float64, p, n)

    for i=1:n
        ind = 0
        for col=1:m
            for row=col+1:m
                ind += 1
                out[ind, i] = xor(X[col, i], X[row, i]) ? -1. : 1.
            end
        end
    end
    out
end

function KLIEP_Hessian(θ, Ψy)
  m, n = size(Ψy)

  w = transpose(Ψy) * θ
  w .= exp.(w)
  w ./= mean(w)

  StatsBase.cov(Ψy, weights(w), 2; corrected=false)
end



function _maxabs(a, b)
    maximum(x -> abs(x[1]-x[2]), zip(a, b))
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
