function SparKLIE_stderr(Ψx, Ψy, θ, ω)
    nx = size(Ψx, 2)
    ny = size(Ψy, 2)
    supp = findall(!iszero, ω)
    S = (cov(transpose(Ψx[supp, :]), corrected=false) ./ nx) .+ (cov(rhat(θ, Ψy) .* transpose(Ψy[supp, :]), corrected=false) ./ ny)
    σ2 = 0.
    for k in 1:length(supp)
        for l in 1:length(supp)
            σ2 += S[k,l] * ω[supp][k] * ω[supp][l]
        end
    end
    sqrt(σ2)
end
