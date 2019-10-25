abstract type KLIEPSolver end

struct CD_KLIEP <: KLIEPSolver end

### interface


KLIEP(Ψx, Ψy, ::CD_KLIEP) = KLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy)

function KLIEP!(x::SparseIterate, Ψx, Ψy)
    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxZero()

    coordinateDescent!(x, f, g)
end

spKLIEP(Ψx, Ψy, λ, ::CD_KLIEP) =
  spKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ)

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

function spKLIEP_refit!(
    x::SparseIterate,
    Ψx::Matrix{Float64},
    Ψy::Matrix{Float64},
    supp::Vector{Int64})

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
    e[row] = -1.0

    x = SparseIterate(m)
    x[row] = 1.
    f = CDQuadraticLoss(H, e)

    # init sigma
    σ = sqrt( dot(x, H * x) )

    for iter=1:10
        coordinateDescent!(x, f, ProxL1(λ0 * σ))
        σnew = sqrt( dot(x, H * x) )

        if abs(σnew - σ) / σ < 1e-3
          break
        end
        σ = σnew
    end

    x
end

function Hinv_row_refit!(
    x::SparseIterate,
    H::Matrix{Float64},
    idx::Int64,
    supp::Vector{Int64})

    m = length(x)
    δ = zeros(m)
    δ[idx] = -1.0

    w = ones(Float64, m) * 1e10
    for k in supp
        w[k] = 0.
    end

    f = CDQuadraticLoss(H, δ)
    g = ProxL1(1., w)
    coordinateDescent!(x, f, g)
end
