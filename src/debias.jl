function debias_KLIEP(idx::Int64, θ, ω, Ψx, Ψy)

    μx = vec(mean(Ψx, dims=2))
    wy = zeros(size(Ψy, 2))
    mul!(wy, transpose(Ψy), θ)
    wy .= exp.(wy)
    wy ./= mean(wy)

    θ1 = θ[idx]
    for k in findall(!iszero, ω)
        θ1 += ω[k] * ( μx[k] - mean( wy .* Ψy[k, :] ) )
    end

    supp = union(idx, findall(!iszero, θ), findall(!iszero, ω))
    spKLIEP_refit!(θ, Ψx, Ψy, supp)

    θ1, θ[idx]
end
