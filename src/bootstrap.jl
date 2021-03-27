struct BootstrapEstimates
    θhat::Vector{Float64}
    θb::Matrix{Float64}     # each column represents a bootstrap estimate
end

function _fill_bΨ!(bΨ, Ψ, b_ind)
    p, n = size(Ψ)
    for col = 1:n
        for row = 1:p
            bΨ[row, col] = Ψ[row, b_ind[col]]
        end
    end
end

function _boot_SparKLIE1(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
    bθ = Matrix{Float64}(undef, length(θ_ind), size(x_ind, 2))
    bΨx = similar(Ψx)
    bΨy = similar(Ψy)
    for b = 1:size(x_ind, 2)
        _fill_bΨ!(bΨx, Ψx, view(x_ind, :, b))
        _fill_bΨ!(bΨy, Ψy, view(y_ind, :, b))
        bμx = vec(mean(bΨx, dims=2))
        for k = 1:length(θ_ind)
            supp = _find_supp(θ_ind[k], θ)
            bθk = KLIEP(bΨx[supp, :], bΨy[supp, :], CD_KLIEP())
            br = rhat(bθk, bΨy[supp, :])
            bθ[k, b] = bθk[end]
            for l in findall(!iszero, Hinv[k])
                bθ[k, b] += Hinv[k][l] * ( bμx[l] - mean( br .* bΨy[l, :] ) )
            end
        end
    end
    bθ
end

function _boot_SparKLIE2(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
    bθ = Matrix{Float64}(undef, length(θ_ind), size(x_ind, 2))
    for k = 1:length(θ_ind)
        supp = _find_supp(θ_ind[k], θ, Hinv[k])
        Ψxk = Ψx[supp, :]
        Ψyk = Ψy[supp, :]
        bΨxk = similar(Ψxk)
        bΨyk = similar(Ψyk)
        for b = 1:size(x_ind, 2)
            _fill_bΨ!(bΨxk, Ψxk, view(x_ind, :, b))
            _fill_bΨ!(bΨyk, Ψyk, view(y_ind, :, b))
            bθk = KLIEP(bΨxk, bΨyk, CD_KLIEP())
            bθ[k, b] = bθk[end]
        end
    end
    bθ
end

function _boot_SparKLIE(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange}; bootSamples, debias)
    if !(size(Ψx, 1) == size(Ψy, 1) == length(θ))
        throw(DimensionMismatch("size(Ψx, 1) = $(size(Ψx, 1)), size(Ψy, 1) = $(size(Ψy, 1)), and length(θ) = $(length(θ)) are not all equal"))
    end
    if length(Hinv) !== length(θ_ind)
        throw(DimensionMismatch("length of Hinv, $(length(Hinv)), does not equal number of parameters, $(length(θ_ind))"))
    end
    nx = size(Ψx, 2)
    ny = size(Ψy, 2)
    x_ind = Matrix{Int16}(undef, nx, bootSamples)
    y_ind = Matrix{Int16}(undef, ny, bootSamples)
    for b = 1:bootSamples
        sample!(1:nx, view(x_ind, :, b))
        sample!(1:ny, view(y_ind, :, b))
    end
    if debias == 1
        θ1 = _debias1(Ψx, Ψy, θ, Hinv, θ_ind)
        bθ1 = _boot_SparKLIE1(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        return BootstrapEstimates(θ1, bθ1)
    elseif debias == 2
        θ2 = _debias2(Ψx, Ψy, θ, Hinv, θ_ind)
        bθ2 = _boot_SparKLIE2(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        return BootstrapEstimates(θ2, bθ2)
    else
        θ1 = _debias1(Ψx, Ψy, θ, Hinv, θ_ind)
        bθ1 = _boot_SparKLIE1(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        θ2 = _debias2(Ψx, Ψy, θ, Hinv, θ_ind)
        bθ2 = _boot_SparKLIE2(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        return BootstrapEstimates(θ1, bθ1), BootstrapEstimates(θ2, bθ2)
    end
end

boot_SparKLIE(Ψx, Ψy, θ, Hinv; bootSamples::Int=300, debias::Int=0) = _boot_SparKLIE(Ψx, Ψy, θ, Hinv, 1:length(θ); bootSamples, debias)
boot_SparKLIE(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange}; bootSamples::Int=300, debias::Int=0) = _boot_SparKLIE(Ψx, Ψy, θ, Hinv, θ_ind; bootSamples, debias)

function boot_quantile(straps::BootstrapEstimates, prob)
    p, bootSamples = size(straps.θb)
    infNormDist = Vector{Float64}(undef, bootSamples)
    for b = 1:bootSamples
        infNormDist[b] = norm_diff(straps.θhat, view(straps.θb, :, b), Inf)
    end
    if count(isnan, infNormDist) == 0
        quantile(infNormDist, prob)
    else
        @warn "NaNs detected in BootstrapEstimates"
        quantile(infNormDist[findall(!isnan, infNormDist)], prob)
    end
end

function boot_quantile_studentized(straps::BootstrapEstimates, prob)
    p, bootSamples = size(straps.θb)
    infNormDist = Vector{Float64}(undef, bootSamples)
    tmp = Vector{Float64}(undef, p)
    w = reshape(std(straps.θb; dims = 2, corrected = false), :)
    for b = 1:bootSamples
        tmp .= (straps.θhat .- straps.θb[:, b]) ./ w
        infNormDist[b] = norm(tmp, Inf)
    end
    if count(isnan, infNormDist) == 0
        quantile(infNormDist, prob), w
    else
        @warn "NaNs detected in BootstrapEstimates"
        quantile(infNormDist[findall(!isnan, infNormDist)], prob), w
    end
end
