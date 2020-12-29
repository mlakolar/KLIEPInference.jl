function KLIEP_stderr(Ψx, Ψy, θ, ω)
    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    supp = findall(!iszero, ω)
    s = length(supp)

    r = zeros(ny)
    mul!(r, transpose(Ψy), θ)
    r .= exp.(r)
    r ./= mean(r)

    Ψyr = Ψy[supp,:]
    for j = 1:ny
        for k = 1:s
            Ψyr[k,j] *= r[j]
        end
    end

    S  = cov( transpose(Ψx[supp,:]), corrected=false) / nx
    S += cov( transpose(Ψyr), corrected=false) / ny

    σ2 = 0.
    for k in 1:s
        for l in 1:s
            σ2 += S[k,l] * ω[supp][k] * ω[supp][l]
        end
    end

    sqrt(σ2)
end
