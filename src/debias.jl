function KLIEP_debias1(Ψx, Ψy, θ, ω, idx::Int64)
    μx = vec(mean(Ψx, dims=2))

    r = zeros(size(Ψy, 2))
    mul!(r, transpose(Ψy), θ)
    r .= exp.(r)
    r ./= mean(r)

    θ1 = θ[idx]
    for k in findall(!iszero, ω)
        θ1 += ω[k] * ( μx[k] - mean( r .* Ψy[k, :] ) )
    end

    θ1
end

function KLIEP_debias1(Ψx, Ψy, θ, ω, idx::Int64, g)
    μx = vec(mean(Ψx, dims=2))

    r = zeros(size(Ψy, 2))
    mul!(r, transpose(Ψy), θ)
    r .= exp.(r)
    r ./= mean(r)

    θ1 = θ[idx]
    for k in findall(!iszero, ω)
        θ1 += ω[k] * ( μx[k] - mean( r .* Ψy[k, :] ) - g[k] )
    end

    θ1
end

function KLIEP_debias2(Ψx, Ψy, θ, ω, idx::Int64)
    supp = sort(union(idx, findall(!iszero, θ), findall(!iszero, ω)))

    θ2 = KLIEP(Ψx[supp,:], Ψy[supp,:], CD_KLIEP())

    θ2[findfirst(x -> x == idx, supp)]
end
