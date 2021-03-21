abstract type KLIEPSolver end

struct CD_KLIEP <: KLIEPSolver end

### interface

KLIEP(Ψx, Ψy, ::CD_KLIEP) = KLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy)

function KLIEP!(x::SparseIterate, Ψx, Ψy)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxZero()
    coordinateDescent!(x, f, g)
end

function spKLIEP!(x::SparseIterate, Ψx, Ψy, λ)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxL1(λ)
    coordinateDescent!(x, f, g)
end

spKLIEP!(x::SparseIterate, f::CDKLIEPLoss, g::ProxL1) = coordinateDescent!(x, f, g)

function _compute_loadings(x, Ψx, Ψy)
    p, nx = size(Ψx)
    ny = size(Ψy, 2)
    Ψyr = transpose( rhat(x, Ψy) .* transpose(Ψy) )
    s = zeros(p)
    for k = 1:p
        s[k] = sqrt( (var(Ψx[k, :]) / nx) + (var(Ψyr[k, :]) / ny) )
    end
    s
end

function spKLIEP_scaled!(x::SparseIterate, Ψx, Ψy, λ)
    f = CDKLIEPLoss(Ψx, Ψy)

    γ = _compute_loadings(x, Ψx, Ψy)
    g = ProxL1(λ, γ)
    coordinateDescent!(x, f, g)

    γ = _compute_loadings(x, Ψx, Ψy)
    g = ProxL1(λ, γ)
    coordinateDescent!(x, f, g, CDOptions(; warmStart=true))
end

spKLIEP(Ψx, Ψy, λ, ::CD_KLIEP; loadings=true) = 
    (loadings ? spKLIEP_scaled!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ) : spKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ))

# Invert row of a Hessian using the scaled lasso
function Hinv_row(H, row, λ0)
    p = size(H, 1)
    e = zeros(p)
    e[row] = -1.

    x = SparseIterate(p)
    x[row] = 1.
    f = CDQuadraticLoss(H, e)

    # init sigma
    σ = sqrt( dot(x, H * x) )
    λ = e .+ 1.

    for iter = 1:10
        coordinateDescent!(x, f, ProxL1(λ0 * σ, λ))
        σnew = sqrt( dot(x, H * x) )
        if abs(σnew - σ) / σ < 1e-3
            break
        end
        σ = σnew
    end

    x
end
