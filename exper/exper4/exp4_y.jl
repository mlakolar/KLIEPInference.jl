# single-edge for base y, a grid of λ1

graphdir = ARGS[1]
graphtype = ARGS[2]
idx = parse(Int, ARGS[3])
rep = parse(Int, ARGS[4])

using KLIEPInference
using JLD
using LinearAlgebra
using SparseArrays

using ProximalBase, CoordinateDescent
include("/home/bkim6/code/refit_to_supp.jl")

file = jldopen("$(graphdir)/$(graphtype)/$(graphtype)_rep_$(rep).jld", "r")

X = read(file, "X")
Y = read(file, "Y")
close(file)

Ψx = Ψising(X .== 1.)
Ψy = Ψising(Y .== 1.)

p, nx = size(Ψx)
ny = size(Ψy, 2)
@show p, nx, ny

λ1 = round.(exp.(range(log(sqrt(4. * log(p) / nx)), stop=log(sqrt(2. * log(p) / nx)), length=5)), digits=3)
λ2 = sqrt(2. * log(p) / ny)

# rows : one-step vs double selection
# cols : λ1
θhat = zeros(2, 5)
σhat = zeros(2, 5)

for t = 1:5
    @show λ1[t]

    θtemp = spKLIEP(Ψx, Ψy, λ1[t], CD_KLIEP())
    θtemp = convert(SparseVector, θtemp)
    supp = union(idx, θtemp.nzind)
    θtemp = SparseIterate(θtemp)
    KLIEP_refit!(θtemp, Ψx, Ψy, supp)

    H = KLIEP_Hessian(θtemp, Ψy)
    ω = Hinv_row(H, idx, λ2)
    ω = convert(SparseVector, ω)
    supp = union(idx, ω.nzind)
    ω = SparseIterate(ω)
    Hinv_row_refit!(ω, H, idx, supp)

    θtemp = convert(SparseVector, θtemp)
    ω = convert(SparseVector, ω)

    θhat[1, t] = KLIEP_debias(idx, θtemp, ω, Ψx, Ψy)
    σhat[1, t] = KLIEP_var(Ψx, Ψy, θtemp, ω)

    supp = union(idx, θtemp.nzind, ω.nzind)
    θtemp = SparseIterate(θtemp)
    KLIEP_refit!(θtemp, Ψx, Ψy, supp)
    ω = SparseIterate(ω)
    Hinv_row_refit!(ω, H, idx, supp)

    θhat[2, t] = θtemp[idx]
    σhat[2, t] = KLIEP_var(Ψx[supp,:], Ψy[supp,:], θtemp[supp], ω[supp])
end

@show θhat

if ~ispath("/home/bkim6/single/$(graphtype)/")
    mkpath("/home/bkim6/single/$(graphtype)/")
end
@save "/home/bkim6/single/$(graphtype)/exp4_y_rep_$(rep).jld" θhat σhat
