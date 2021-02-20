struct CD_SymKLIEP <: KLIEPSolver end

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
    mul!(f.rx, transpose(f.Ψx), x)
    f.rx .*= -1.

    mul!(f.ry, transpose(f.Ψy), x)

    nothing
end

function CoordinateDescent.gradient(
    f::CDSymKLIEPLoss,
    x::SparseIterate,
    j::Int64)

    μx, Ψx, rx = f.μx, f.Ψx, f.rx
    μy, Ψy, ry = f.μy, f.Ψy, f.ry

    -μx[j] + μy[j] - mean( exp.(rx) .* Ψx[j, :] ) / mean(exp, rx) + mean( exp.(ry) .* Ψy[j, :] ) / mean(exp, ry)
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
    mean_exp_rx_Ψx2 = mean( exp.(rx) .* Ψx[j, :].^2 )

    mean_exp_ry = mean(exp, ry)
    mean_exp_ry_Ψy = mean( exp.(ry) .* Ψy[j, :] )
    mean_exp_ry_Ψy2 = mean( exp.(ry) .* Ψy[j, :].^2 )

    @inbounds grad = -μx[j] + μy[j] - mean_exp_rx_Ψx / mean_exp_rx + mean_exp_ry_Ψy / mean_exp_ry
    H = mean_exp_rx_Ψx2 / mean_exp_rx - mean_exp_rx_Ψx^2 / mean_exp_rx^2 + mean_exp_ry_Ψy2 / mean_exp_ry - mean_exp_ry_Ψy^2 / mean_exp_ry^2

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
    m, ny = size(Ψy)

    wx = -transpose(Ψx) * θ
    wx .= exp.(wx)
    wx ./= mean(wx)

    wy = transpose(Ψy) * θ
    wy .= exp.(wy)
    wy ./= mean(wy)

    StatsBase.cov(Ψx, weights(wx), 2; corrected=false) + StatsBase.cov(Ψy, weights(wy), 2; corrected=false)
end

####################################
#
# symmetric KLIEP loss
#
####################################

SymKLIEP(Ψx, Ψy, ::CD_SymKLIEP) = SymKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy)

function SymKLIEP!(x::SparseIterate, Ψx, Ψy)
    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxZero()

    coordinateDescent!(x, f, g)
end

spSymKLIEP(Ψx, Ψy, λ, ::CD_SymKLIEP) = spSymKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ)

function spSymKLIEP!(x::SparseIterate, Ψx, Ψy, λ)
    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxL1(λ)

    coordinateDescent!(x, f, g)
end

function spSymKLIEP_refit!(x::SparseIterate, Ψx, Ψy)
    m = length(x)
    ω = ones(Float64, m) * 1e10

    for i=1:m
        if !iszero(x[i])
            ω[i] = 0.
        end
    end

    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxL1(1., ω)
    coordinateDescent!(x, f, g)
end

function spSymKLIEP_refit!(
    x::SparseIterate,
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},
    supp::Vector{Int64})

    w = ones(Float64, length(x)) * 1e10
    for k in supp
        w[k] = 0.
    end

    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxL1(1., w)
    coordinateDescent!(x, f, g)
end

function SymKLIEP_debias(
    ind::Int64,
    θ,
    ω::SparseVector,
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64})

    μx = vec(mean(Ψx, dims=2))
    μy = vec(mean(Ψy, dims=2))

    wx = zeros(size(Ψx,2))
    mul!(wx, transpose(Ψx), θ)
    wx .*= -1.
    wx .= exp.(wx)
    wx ./= mean(wx)

    wy = zeros(size(Ψy,2))
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    θ1 = θ[ind]
    for k in ω.nzind
        θ1 += ω[k] * ( μx[k] - μy[k] + mean( wx .* Ψx[k, :] ) - mean( wy .* Ψy[k, :] ) )
    end
    return θ1
end

# variance: sparse ω
function SymKLIEP_var(
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},
    θ,
    ω::SparseVector)

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    supp = ω.nzind
    s = nnz(ω)

    # compute scalars
    wx = -transpose(Ψx) * θ
    wx .= exp.(wx)
    wx ./= mean(wx)
    wx .+= 1.

    wy = transpose(Ψy) * θ
    wy .= exp.(wy)
    wy ./= mean(wy)
    wy .+= 1.

    # compute scaled data
    Ψxw = Ψx[supp, :]
    for i = 1:nx
        for k = 1:s
            Ψxw[k, i] *= wy[i]
        end
    end

    Ψyw = Ψy[supp, :]
    for j = 1:ny
        for k = 1:s
            Ψyw[k, j] *= wy[j]
        end
    end

    # compute sample covariances
    S  = cov( transpose(Ψxw), corrected=false ) / nx
    S += cov( transpose(Ψyw), corrected=false ) / ny

    # combine
    σhat2 = 0.
    for k in 1:s
        for l in 1:s
            σhat2 += S[k,l] * ω.nzval[k] * ω.nzval[l]
        end
    end

    return σhat2
end

# variance: dense ω
function SymKLIEP_var(
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},
    θ,
    ω::Vector{Float64})

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    # compute scalars
    wx = -transpose(Ψx) * θ
    wx .= exp.(wx)
    wx ./= mean(wx)
    wx .+= 1.

    wy = transpose(Ψy) * θ
    wy .= exp.(wy)
    wy ./= mean(wy)
    wy .+= 1.

    # compute scaled data
    Ψxw = zeros(p,nx)
    for i = 1:nx
        for k = 1:p
            Ψxw[k, i] = Ψx[k, i] * wx[i]
        end
    end

    Ψyw = zeros(p,ny)
    for j = 1:ny
        for k = 1:p
            Ψyw[k, j] = Ψy[k, j] * wy[j]
        end
    end

    # compute sample covariances
    S  = cov( transpose(Ψxw), corrected=false ) / nx
    S += cov( transpose(Ψyw), corrected=false ) / ny

    # combine
    σhat2 = 0.
    for k in 1:p
        for l in 1:p
            σhat2 += S[k,l] * ω[k] * ω[l]
        end
    end

    return σhat2
end
