# balanced sample sizes, single-edge for base y, a grid of λ1

graphdir = ARGS[1]
graphtype = ARGS[2]
idx = parse(Int, ARGS[3])
nx = parse(Int, ARGS[4])
ny = parse(Int, ARGS[5])
rep = parse(Int, ARGS[6])

using KLIEPInference
using ProximalBase, CoordinateDescent
using JLD
using LinearAlgebra
using SparseArrays

include("/home/bkim6/code/refit_to_supp.jl")

file = jldopen("$(graphdir)/$(graphtype)/$(graphtype)_rep_$(rep).jld", "r")
X = read(file, "X")
Y = read(file, "Y")
close(file)

X = X[:, 1:nx]
Y = Y[:, 1:ny]

Ψx = Ψising(X .== 1.)
Ψy = Ψising(Y .== 1.)

p = size(Ψx, 1)

@show p, nx, ny

λ1 = round.(exp.(range(log(sqrt(16. * log(p) / nx)), stop=log(sqrt(2. * log(p) / nx)), length=5)), digits=3)
λ2 = sqrt(2. * log(p) / ny)

# rows : one-step vs double selection
# cols : λ1
θhat = zeros(2, 5)
σhat = zeros(2, 5)

for t = 1:5
    @show λ1[t]

    θ = spKLIEP(Ψx, Ψy, λ1[t], CD_KLIEP())
    θ = convert(SparseVector, θ)
    supp = union(idx, θ.nzind)
    θ = SparseIterate(θ)
    spKLIEP_refit!(θ, Ψx, Ψy, supp)

    H = KLIEP_Hessian(θ, Ψy)
    ω = Hinv_row(H, idx, λ2)
    ω = convert(SparseVector, ω)
    supp = union(idx, ω.nzind)
    ω = SparseIterate(ω)
    Hinv_row_refit!(ω, H, idx, supp)

    θ = convert(SparseVector, θ)
    ω = convert(SparseVector, ω)

    θhat[1, t] = KLIEP_debias(idx, θ, ω, Ψx, Ψy)
    σhat[1, t] = KLIEP_var(Ψx, Ψy, θ, ω)

    supp = union(idx, θ.nzind, ω.nzind)
    θ = SparseIterate(θ)
    spKLIEP_refit!(θ, Ψx, Ψy, supp)
    ω = SparseIterate(ω)
    Hinv_row_refit!(ω, H, idx, supp)

    θhat[2, t] = θ[idx]
    σhat[2, t] = KLIEP_var(Ψx[supp,:], Ψy[supp,:], θ[supp], ω[supp])
end

@show θhat

if ~ispath("/home/bkim6/single/$(graphtype)/")
    mkpath("/home/bkim6/single/$(graphtype)/")
end
@save "/home/bkim6/single/$(graphtype)/exp4_y_rep_$(rep).jld" θhat σhat
