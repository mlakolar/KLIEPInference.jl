function KLIEP_debias(
    ind::Int64,
    θ,
    ω::SparseVector,
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},)

    μx = vec(mean(Ψx, dims=2))
    wy = zeros(size(Ψy,2))
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    θ1 = θ[ind]
    for k in ω.nzind
        θ1 += ω[k] * ( μx[k] - mean( wy .* Ψy[k, :] ) )
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
    for j = 1:ny
        for k = 1:s
            Ψysub[k, j] *= wy[j]
        end
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
    for j = 1:ny
        for k = 1:p
            Ψyw[k, j] = wy[j] * Ψy[k, j]
        end
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
