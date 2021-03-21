function _debias1(Ψx, Ψy, θ, ω, θ_ind::Int; refit=true)
    μx = vec(mean(Ψx, dims=2))
    if refit == false
        r = rhat(θ, Ψy)
        θ1 = θ[θ_ind]
        for l in findall(!iszero, ω)
            θ1 += ω[l] * ( μx[l] - mean( r .* Ψy[l, :] ) )
        end
        return θ1
    else
        supp = _find_supp(θ_ind, θ)
        θk = KLIEP(Ψx[supp, :], Ψy[supp, :], CD_KLIEP())

        r = rhat(θk, Ψy[supp, :])
        θ1 = θk[end]
        for l in findall(!iszero, ω)
            θ1 += ω[l] * ( μx[l] - mean( r .* Ψy[l, :] ) )
        end
        return θ1
    end
end

function _debias1(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange}; refit=true)
    θ1 = Vector{Float64}(undef, length(θ_ind))
    μx = vec(mean(Ψx, dims=2))
    if refit == false
        for k in 1:length(θ_ind)
            r = rhat(θ, Ψy)
            θ1[k] = θk[θ_ind[k]]
            for l in findall(!iszero, Hinv[k])
                θ1[k] += Hinv[k][l] * ( μx[l] - mean( r .* Ψy[l, :] ) )
            end
        end
        return θ1
    else
        for k in 1:length(θ_ind)
            supp = _find_supp(θ_ind[k], θ)
            θk = KLIEP(Ψx[supp, :], Ψy[supp, :], CD_KLIEP())

            r = rhat(θk, Ψy[supp, :])
            θ1[k] = θk[end]
            for l in findall(!iszero, Hinv[k])
                θ1[k] += Hinv[k][l] * ( μx[l] - mean( r .* Ψy[l, :] ) )
            end
        end
        return θ1
    end
end

function _debias2(Ψx, Ψy, θ, ω, θ_ind::Int)
    supp = _find_supp(θ_ind, θ, ω)
    θk = KLIEP(Ψx[supp,:], Ψy[supp,:], CD_KLIEP())
    θk[end]
end

function _debias2(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange})
    θ2 = Vector{Float64}(undef, length(θ_ind))
    for k in 1:length(θ_ind)
        supp = _find_supp(θ_ind[k], θ, Hinv[k])
        θk = KLIEP(Ψx[supp, :], Ψy[supp, :], CD_KLIEP())
        θ2[k] = θk[end]
    end
    θ2
end
