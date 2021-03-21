####################################
#
# symmetric KLIEP loss
#
####################################

struct CDSymKLIEPLoss <: CoordinateDifferentiableFunction
    μx::Vector{Float64}   # mean(Ψx, dims = 2)
    μy::Vector{Float64}   # mean(Ψy, dims = 2)
    Ψx::Matrix{Float64}
    Ψy::Matrix{Float64}
    rx::Vector{Float64}    # nx dim vector --- stores -θ'Ψx(i)
    ry::Vector{Float64}    # ny dim vector --- stores  θ'Ψy(j)
    p::Int64
end

function CDSymKLIEPLoss(Ψx::Matrix{Float64}, Ψy::Matrix{Float64})
    (p = size(Ψx, 1)) == size(Ψy, 1) || throw(DimensionMismatch())

    CDSymKLIEPLoss(vec(mean(Ψx, dims=2)), vec(mean(Ψy, dims=2)), Ψx, Ψy, zeros(size(Ψx, 2)), zeros(size(Ψy, 2)), p)
end

CoordinateDescent.numCoordinates(f::CDSymKLIEPLoss) = f.p

function CoordinateDescent.initialize!(f::CDSymKLIEPLoss, x::SparseIterate)
    mul!(f.rx, transpose(f.Ψx), x)
    f.rx .*= -1.

    mul!(f.ry, transpose(f.Ψy), x)

    nothing
end

function CoordinateDescent.gradient(
    f::CDSymKLIEPLoss,
    x::SparseIterate,
    j::Int64)

    μx, Ψx, rx = f.μx, f.Ψx, f.rx
    μy, Ψy, ry = f.μy, f.Ψy, f.ry

    -μx[j] + μy[j] - mean( exp.(rx) .* Ψx[j, :] ) / mean(exp, rx) + mean( exp.(ry) .* Ψy[j, :] ) / mean(exp, ry)
end

function CoordinateDescent.descendCoordinate!(
    f::CDSymKLIEPLoss,
    g::Union{ProxL1, ProxZero},
    x::SparseIterate,
    j::Int64)

    μx, Ψx, rx = f.μx, f.Ψx, f.rx
    μy, Ψy, ry = f.μy, f.Ψy, f.ry

    mean_exp_rx = mean(exp, rx)
    mean_exp_rx_Ψx = mean( exp.(rx) .* Ψx[j, :] )
    mean_exp_rx_Ψx2 = mean( exp.(rx) .* Ψx[j, :].^2 )

    mean_exp_ry = mean(exp, ry)
    mean_exp_ry_Ψy = mean( exp.(ry) .* Ψy[j, :] )
    mean_exp_ry_Ψy2 = mean( exp.(ry) .* Ψy[j, :].^2 )

    @inbounds grad = -μx[j] + μy[j] - mean_exp_rx_Ψx / mean_exp_rx + mean_exp_ry_Ψy / mean_exp_ry
    H = mean_exp_rx_Ψx2 / mean_exp_rx - mean_exp_rx_Ψx^2 / mean_exp_rx^2 + mean_exp_ry_Ψy2 / mean_exp_ry - mean_exp_ry_Ψy^2 / mean_exp_ry^2

    @inbounds oldVal = x[j]
    @inbounds x[j] -= grad / H
    newVal = cdprox!(g, x, j, 1. / H)
    h = newVal - oldVal

    # update internals
    for i=1:length(rx)
        @inbounds rx[i] -= h*Ψx[j, i]
    end
    for i=1:length(ry)
        @inbounds ry[i] += h*Ψy[j, i]
    end
    h
end


# symmetric KLIEP solver
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


# Computes the Hessian of the SymKLIEP loss
SymKLIEP_Hessian(θ, Ψx, Ψy) = StatsBase.cov(Ψx, weights(KLIEPInference.rhat(-θ, Ψx)), 2; corrected=false) + StatsBase.cov(Ψy, weights(KLIEPInference.rhat(θ, Ψy)), 2; corrected=false)


# debiasing
function SymKLIEP_debias1(Ψx, Ψy, θ, ω, θ_ind::Int)
    μx = vec(mean(Ψx, dims=2))
    μy = vec(mean(Ψy, dims=2))

    supp = KLIEPInference._find_supp(θ_ind, θ)
    θk = SymKLIEP(Ψx[supp, :], Ψy[supp, :], CD_SymKLIEP())

    rx = KLIEPInference.rhat(-θk, Ψx[supp, :])
    ry = KLIEPInference.rhat( θk, Ψy[supp, :])
    θ1 = θk[end]
    for l in findall(!iszero, ω)
        θ1 += ω[l] * ( μx[l] - μy[l] + mean( rx .* Ψx[l, :] ) - mean( ry .* Ψy[l, :] ) )
    end
    θ1
end

function SymKLIEP_debias2(Ψx, Ψy, θ, ω, θ_ind::Int)
    supp = KLIEPInference._find_supp(θ_ind, θ, ω)
    θk = SymKLIEP(Ψx[supp,:], Ψy[supp,:], CD_SymKLIEP())
    θk[end]
end


# standard error
function SymKLIEP_stderr(Ψx, Ψy, θ, ω)
    nx = size(Ψx, 2)
    ny = size(Ψy, 2)
    supp = findall(!iszero, ω)
    S = (cov(transpose(Ψx[supp, :]), corrected=false) ./ nx) .+ (cov(transpose(Ψy[supp, :]), corrected=false) ./ ny) .+ (cov(KLIEPInference.rhat(-θ, Ψx) .* transpose(Ψx[supp, :]), corrected=false) ./ nx) .+ (cov(KLIEPInference.rhat(θ, Ψy) .* transpose(Ψy[supp, :]), corrected=false) ./ ny)
    σ2 = 0.
    for k in 1:length(supp)
        for l in 1:length(supp)
            σ2 += S[k,l] * ω[supp][k] * ω[supp][l]
        end
    end
    sqrt(σ2)
end
