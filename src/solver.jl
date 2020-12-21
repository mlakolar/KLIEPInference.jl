abstract type KLIEPSolver end

struct CD_KLIEP <: KLIEPSolver end

### interface


KLIEP(Ψx, Ψy, ::CD_KLIEP) = KLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy)

function KLIEP!(x::SparseIterate, Ψx, Ψy)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxZero()
    coordinateDescent!(x, f, g)
end


function spKLIEP(Ψx, Ψy, λ, ::CD_KLIEP; loadings=true)
    if loadings
        spKLIEP_scaled!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ)
    else
        spKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ)
    end
end

function _compute_loadings(x, Ψx, Ψy)
    p, nx = size(Ψx)
    ny = size(Ψy, 2)

    wy = zeros(ny)
    mul!(wy, transpose(Ψy), x)
    wy .= exp.(wy)
    wy ./= mean(wy)

    Ψyw = zeros(p, ny)
    for j = 1:ny
        for k = 1:p
            Ψyw[k,j] = wy[j] * Ψy[k,j]
        end
    end

    s = zeros(p)
    for k = 1:p
        s[k] = sqrt((var(Ψx[k,:]) / nx) + (var(Ψyw[k,:]) / ny))
    end

    return s
end

function spKLIEP_scaled!(x::SparseIterate, Ψx, Ψy, λ)
    f = CDKLIEPLoss(Ψx, Ψy)

    γ = _compute_loadings(x, Ψx, Ψy)
    g = ProxL1(λ, γ)
    coordinateDescent!(x, f, g)

    γ = _compute_loadings(x, Ψx, Ψy)
    g = ProxL1(λ, γ)
    coordinateDescent!(x, f, g, CDOptions(; warmStart=true))

    return x
end

function spKLIEP!(x::SparseIterate, Ψx, Ψy, λ)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxL1(λ)
    coordinateDescent!(x, f, g)
end

spKLIEP!(x::SparseIterate, f::CDKLIEPLoss, g::ProxL1) = coordinateDescent!(x, f, g)

function spKLIEP_refit!(x::SparseIterate, Ψx, Ψy)
    m = length(x)
    ω = ones(Float64, m) * 1e10

    for i=1:m
        if !iszero(x[i])
            ω[i] = 0.
        end
    end

    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxL1(1., ω)
    coordinateDescent!(x, f, g)
end

function spKLIEP_refit!(x::SparseIterate, Ψx, Ψy, supp::Vector{Int64})
    w = ones(Float64, length(x)) * 1e10

    for k in supp
        w[k] = 0.
    end

    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxL1(1., w)
    coordinateDescent!(x, f, g)
end


function Hinv_row(H, row, λ0)
    m = size(H, 1)
    e = zeros(m)
    e[row] = -1.

    x = SparseIterate(m)
    x[row] = 1.
    f = CDQuadraticLoss(H, e)

    # init sigma
    σ = sqrt( dot(x, H * x) )
    λ = e .+ 1.

    for iter=1:10
        coordinateDescent!(x, f, ProxL1(λ0 * σ, λ))
        σnew = sqrt( dot(x, H * x) )

        if abs(σnew - σ) / σ < 1e-3
            break
        end
        σ = σnew
    end

    x
end

function Hinv_row_refit!(x::SparseIterate, H, row, supp::Vector{Int64})
    m = length(x)
    δ = zeros(m)
    δ[row] = -1.0

    w = ones(Float64, m) * 1e10
    for k in supp
        w[k] = 0.
    end

    f = CDQuadraticLoss(H, δ)
    g = ProxL1(1., w)
    coordinateDescent!(x, f, g)
end
