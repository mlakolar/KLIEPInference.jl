function KLIEP_debias(
    ind::Int64,
    θ,
    ω::SparseVector,
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},)

    μx = vec(mean(Ψx, dims=2))
    ry = zeros(size(Ψy,2))
    mul!(ry, transpose(Ψy), θ)

    θ1 = θ[ind]
    for k in ω.nzind
        θ1 += ω[k] * (μx[k] - mean( exp.(ry) .* Ψy[k, :] ) / mean(exp, ry))
    end
    return θ1
end

function SymKLIEP_debias(
    ind::Int64,
    θ,
    ω::SparseVector,
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},)

    μx = vec(mean(Ψx, dims=2))
    μy = vec(mean(Ψy, dims=2))

    rx = zeros(size(Ψx,2))
    mul!(rx, transpose(Ψx), θ)
    rx .*= -1.

    ry = zeros(size(Ψy,2))
    mul!(ry, transpose(Ψy), θ)

    θ1 = θ[ind]
    for k in ω.nzind
        θ1 += ω[k] * (μx[k] - μy[k] + mean( exp.(rx) .* Ψx[k, :] ) / mean(exp, rx) - mean( exp.(ry) .* Ψy[k, :] ) / mean(exp, ry))
    end
    return θ1
end

# variance: sparse ω
function KLIEP_var(
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},
    θ,
    ω::SparseVector)

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    supp = ω.nzind
    s = nnz(ω)

    # compute scalars
    wy = zeros(ny)
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    # compute scaled data
    Ψxsub = Ψx[supp, :]
    Ψysub = Ψy[supp, :]
    for k = 1:s
        Ψysub[k, :] *= wy[k]
    end

    # compute sample covariances
    S  = cov( transpose(Ψxsub), corrected=false) / nx
    S += cov( transpose(Ψysub), corrected=false) / ny

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
function KLIEP_var(
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},
    θ,
    ω::Vector{Float64})

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    # compute scalars
    wy = zeros(ny)
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    # compute scaled data
    Ψyw = zeros(p,ny)
    for k = 1:p
        Ψyw[k, :] = wy[k] * Ψy[k, :]
    end

    # compute sample covariances
    S  = cov( transpose(Ψx) , corrected=false) / nx
    S += cov( transpose(Ψyw), corrected=false) / ny

    # combine
    σhat2 = 0.
    for k in 1:p
        for l in 1:p
            σhat2 += S[k,l] * ω[k] * ω[l]
        end
    end

    return σhat2
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
    for k = 1:s
        Ψxw[k, :] *= wy[k]
    end

    Ψyw = Ψy[supp, :]
    for k = 1:s
        Ψyw[k, :] *= wy[k]
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
    for k = 1:p
        Ψxw[k, :] = Ψx[k, :] * wx[k]
    end

    Ψyw = zeros(p,ny)
    for k = 1:p
        Ψyw[k, :] = Ψy[k, :] * wy[k]
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
