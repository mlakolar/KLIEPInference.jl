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

function _find_supp3(θ, ω, k)
    supp1 = findall(!iszero, θ)
    supp2 = findall(!iszero, ω)
    supp3 = union(supp1, supp2)
    pos_k = findfirst(isequal(k), supp3)
    if pos_k == nothing
        push!(supp3, k)
    else
        supp3[pos_k], supp3[end] = supp3[end], supp3[pos_k]
    end
    supp3
end

function _boot_SparKLIE1(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
    bΨx = similar(Ψx)
    bΨy = similar(Ψy)
    θb = Matrix{Float64}(undef, length(θ_ind), size(x_ind, 2))
    for b = 1:size(x_ind, 2)
        _fill_bΨ!(bΨx, Ψx, view(x_ind, :, b))
        _fill_bΨ!(bΨy, Ψy, view(y_ind, :, b))
        bμx = vec(mean(bΨx, dims=2))
        for k = 1:length(θ_ind)
            supp3 = _find_supp3(θ, [], θ_ind[k])
            bθk = KLIEP(bΨx[supp3, :], bΨy[supp3, :], CD_KLIEP())
            br = rhat(bθk, bΨy[supp3, :])
            θb[k, b] = bθk[end]
            for l in findall(!iszero, Hinv[k])
                θb[k, b] += Hinv[k][l] * ( bμx[l] - mean( br .* bΨy[l, :] ) )
            end
        end
    end
    θb
end

function _boot_SparKLIE2(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
    θb = Matrix{Float64}(undef, length(θ_ind), size(x_ind, 2))
    for k = 1:length(θ_ind)
        supp3 = _find_supp3(θ, Hinv[k], θ_ind[k])
        Ψxk = Ψx[supp3, :]
        Ψyk = Ψy[supp3, :]
        bΨxk = similar(Ψxk)
        bΨyk = similar(Ψyk)
        for b = 1:size(x_ind, 2)
            _fill_bΨ!(bΨxk, Ψxk, view(x_ind, :, b))
            _fill_bΨ!(bΨyk, Ψyk, view(y_ind, :, b))
            θk = KLIEP(bΨxk, bΨyk, CD_KLIEP())
            θb[k, b] = θk[end]
        end
    end
    θb
end

function _boot_SparKLIE(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange}; bootSamples, debias)
    if !(size(Ψx, 1) === size(Ψy, 1) === length(θ))
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
    if debias === 1
        θhat = _debias1(Ψx, Ψy, θ, Hinv, θ_ind)
        θb = _boot_SparKLIE1(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        return BootstrapEstimates(θhat, θb)
    elseif debias === 2
        θhat = _debias2(Ψx, Ψy, θ, Hinv, θ_ind)
        θb = _boot_SparKLIE2(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        return BootstrapEstimates(θhat, θb)
    else
        θhat1 = _debias1(Ψx, Ψy, θ, Hinv, θ_ind)
        θhat2 = _debias2(Ψx, Ψy, θ, Hinv, θ_ind)
        θb1 = _boot_SparKLIE1(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        θb2 = _boot_SparKLIE2(Ψx, Ψy, θ, Hinv, θ_ind, x_ind, y_ind)
        return BootstrapEstimates(θhat1, θb1), BootstrapEstimates(θhat2, θb2)
    end
end

boot_SparKLIE(Ψx, Ψy, θ, Hinv, ::Nothing; bootSamples::Int=300, debias::Int=0) =
    _boot_SparKLIE(Ψx, Ψy, θ, Hinv, 1:length(θ); bootSamples, debias)
boot_SparKLIE(Ψx, Ψy, θ, Hinv, θ_ind::Union{Vector{Int},UnitRange}; bootSamples::Int=300, debias::Int=0) =
    _boot_SparKLIE(Ψx, Ψy, θ, Hinv, θ_ind; bootSamples, debias)

function simulCI(straps::BootstrapEstimates, α::Float64=0.05)
    m, bootSamples = size(straps.θb)

    infNormDist = Vector{Float64}(undef, bootSamples)
    CI = Matrix{Float64}(undef, m, 2)

    for b=1:bootSamples
        infNormDist[b] = norm_diff(straps.θhat, view(straps.θb, :, b), Inf)
    end
    x = quantile!(infNormDist, 1 - α)
    CI[:, 1] .= straps.θhat .- x
    CI[:, 2] .= straps.θhat .+ x

    CI
end

function simulCIstudentized(straps::BootstrapEstimates, α::Float64=0.05)
    m, bootSamples = size(straps.θb)

    infNormDist = Vector{Float64}(undef, bootSamples)
    w = reshape(std(straps.θb; dims = 2, corrected = false), :)

    CI = Matrix{Float64}(undef, m, 2)
    tmp = Vector{Float64}(undef, m)

    for b=1:bootSamples
        tmp .= (straps.θhat .- straps.θb[:, b]) ./ w
        infNormDist[b] = norm(tmp, Inf)
    end
    x = quantile!(infNormDist, 1 - α)
    @. CI[:, 1] = straps.θhat - x * w
    @. CI[:, 2] = straps.θhat + x * w

    CI
end
