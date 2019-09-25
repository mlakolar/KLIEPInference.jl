# Question: Should CD_SymKLIEP be a different type?

struct CD_SymKLIEP <: KLIEPSolver end

SymKLIEP(Ψx, Ψy, ::CD_SymKLIEP) = SymKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy)

function SymKLIEP!(x::SparseIterate, Ψx, Ψy)
    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxZero()

    coordinateDescent!(x, f, g)
end

spSymKLIEP(Ψx, Ψy, λ, ::CD_SymKLIEP) = spSymKLIEP!(SparseIterate(size(Ψx, 1)), Ψx, Ψy, λ)

function spSymKLIEP!(x::SparseIterate, Ψx, Ψy, λ)
    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxL1(λ)

    coordinateDescent!(x, f, g)
end

function spSymKLIEP_refit!(x::SparseIterate, Ψx, Ψy)
    m = length(x)
    ω = ones(Float64, m) * 1e10

    for i=1:m
        if !iszero(x[i])
            ω[i] = 0.
        end
    end

    f = CDSymKLIEPLoss(Ψx, Ψy)
    g = ProxL1(1., ω)
    coordinateDescent!(x, f, g)
end
