function stderr_KLIEP(θ, ω, Ψx, Ψy)

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    wy = zeros(ny)
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    Ψyw = zeros(p, ny)
    for j = 1:ny
        for k = 1:p
            Ψyw[k, j] = wy[j] * Ψy[k, j]
        end
    end

    S  = cov( transpose(Ψx) , corrected=false) / nx
    S += cov( transpose(Ψyw), corrected=false) / ny

    σ2hat = 0.
    for k in 1:p
        for l in 1:p
            σ2hat += S[k,l] * ω[k] * ω[l]
        end
    end

    sqrt(σ2hat)
end

function stderr_KLIEP(θ, ω::SparseIterate, Ψx, Ψy)

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    supp = findall(!iszero, ω)
    s = length(supp)

    wy = zeros(ny)
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    Ψxsub = Ψx[supp, :]
    Ψysub = Ψy[supp, :]
    for j = 1:ny
        for k = 1:s
            Ψysub[k, j] *= wy[j]
        end
    end

    S  = cov( transpose(Ψxsub), corrected=false) / nx
    S += cov( transpose(Ψysub), corrected=false) / ny

    σ2hat = 0.
    for k in 1:s
        for l in 1:s
            σ2hat += S[k,l] * ω[supp][k] * ω[supp][l]
        end
    end

    sqrt(σ2hat)
end
