####################################
#
# symmetric KLIEP loss
#
####################################
struct CDSymKLIEPLoss <: CoordinateDifferentiableFunction
    μx::Vector{Float64}   # mean(Ψx, dims = 2)
    μy::Vector{Float64}   # mean(Ψy, dims = 2)
    Ψx::Matrix{Float64}
    Ψy::Matrix{Float64}
    rx::Vector{Float64}    # nx dim vector --- stores -θ'Ψx(i)
    ry::Vector{Float64}    # ny dim vector --- stores  θ'Ψy(j)
    p::Int64
end

function CDSymKLIEPLoss(Ψx::Matrix{Float64}, Ψy::Matrix{Float64})
    (p = size(Ψx, 1)) == size(Ψy, 1) || throw(DimensionMismatch())

    CDSymKLIEPLoss(vec(mean(Ψx, dims=2)), vec(mean(Ψy, dims=2)), Ψx, Ψy, zeros(size(Ψx, 2)), zeros(size(Ψy, 2)), p)
end

CoordinateDescent.numCoordinates(f::CDSymKLIEPLoss) = f.p

function CoordinateDescent.initialize!(f::CDSymKLIEPLoss, x::SparseIterate)
    mul!(f.rx, transpose(f.Ψx),-x)
    mul!(f.ry, transpose(f.Ψy), x)

    nothing
end

function CoordinateDescent.gradient(
  f::CDSymKLIEPLoss,
  x::SparseIterate,
  j::Int64)

  μx, Ψx, rx = f.μx, f.Ψx, f.rx
  μy, Ψy, ry = f.μy, f.Ψy, f.ry

  -μx[j] + μy[j] + mean( exp.(rx) .* Ψx[j, :] ) / mean(exp, rx) + mean( exp.(ry) .* Ψy[j, :] ) / mean(exp, ry)
end

function CoordinateDescent.descendCoordinate!(
    f::CDSymKLIEPLoss,
    g::Union{ProxL1, ProxZero},
    x::SparseIterate,
    j::Int64)

    μx, Ψx, rx = f.μx, f.Ψx, f.rx
    μy, Ψy, ry = f.μy, f.Ψy, f.ry

    mean_exp_rx = mean(exp, rx)
    mean_exp_rx_Ψx = mean( exp.(rx) .* Ψx[j, :] )
    mean_exp_rx_Ψx2 = mean( exp.(rx) .* Ψx[j, :].^2. )

    mean_exp_ry = mean(exp, ry)
    mean_exp_ry_Ψy = mean( exp.(ry) .* Ψy[j, :] )
    mean_exp_ry_Ψy2 = mean( exp.(ry) .* Ψy[j, :].^2. )

    @inbounds grad = -μx[j] + μy[j] + mean_exp_rx_Ψx / mean_exp_rx + mean_exp_ry_Ψy / mean_exp_ry
    H = mean_exp_rx_Ψx2 / mean_exp_rx - mean_exp_rx_Ψx^2 / mean_exp_rx^2 + mean_exp_ry_Ψy2 / mean_exp_ry - mean_exp_ry_Ψy^2 / mean_exp_ry^2.

    @inbounds oldVal = x[j]
    @inbounds x[j] -= grad / H
    newVal = cdprox!(g, x, j, 1. / H)
    h = newVal - oldVal

    # update internals
    for i=1:length(rx)
        @inbounds rx[i] -= h*Ψx[j, i]
    end
    for i=1:length(ry)
        @inbounds ry[i] += h*Ψy[j, i]
    end
    h
end


function SymKLIEP_Hessian(θ, Ψx, Ψy)
    m, n = size(Ψy)

    wx = transpose(Ψx) *-θ
    wx .= exp.(wx)
    wx ./= mean(wx)

    wy = transpose(Ψy) * θ
    wy .= exp.(wy)
    wy ./= mean(wy)

    StatsBase.cov(Ψx, weights(wx), 2; corrected=false) + StatsBase.cov(Ψy, weights(wy), 2; corrected=false)
end
