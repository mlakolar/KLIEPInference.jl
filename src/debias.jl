function _debias1(Ψx, Ψy, θ, ω, θ_ind::Int)
    μx = vec(mean(Ψx, dims=2))
    supp3 = _findSupp3(θ, [], θ_ind)
    θk = KLIEP(Ψx[supp3, :], Ψy[supp3, :], CD_KLIEP())
    r = rhat(θk, Ψy[supp3, :])
    θ1 = θk[end]
    for l in findall(!iszero, ω)
        θ1 += ω[l] * ( μx[l] - mean( r .* Ψy[l, :] ) )
    end
    θ1
end

function _debias1(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange})
    θ1 = Vector{Float64}(undef, length(θ_ind))
    μx = vec(mean(Ψx, dims=2))
    for k in 1:length(θ_ind)
        supp3 = _findSupp3(θ, [], θ_ind[k])
        θk = KLIEP(Ψx[supp3, :], Ψy[supp3, :], CD_KLIEP())
        r = rhat(θk, Ψy[supp3, :])
        θ1[k] = θk[end]
        for l in findall(!iszero, Hinv[k])
            θ1[k] += Hinv[k][l] * ( μx[l] - mean( r .* Ψy[l, :] ) )
        end
    end
    θ1
end

function _debias2(Ψx, Ψy, θ, ω, θ_ind::Int)
    supp3 = _findSupp3(θ, ω, θ_ind)
    θk = KLIEP(Ψx[supp3,:], Ψy[supp3,:], CD_KLIEP())
    θk[end]
end

function _debias2(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange})
    θ2 = Vector{Float64}(undef, length(θ_ind))
    for k in 1:length(θ_ind)
        supp3 = _findSupp3(θ, Hinv[k], θ_ind[k])
        θk = KLIEP(Ψx[supp3, :], Ψy[supp3, :], CD_KLIEP())
        θ2[k] = θk[end]
    end
    θ2
end
