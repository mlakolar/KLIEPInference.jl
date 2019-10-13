# single-edge for base y, a grid of λ1

graphtype = ARGS[1]
idx = parse(Int, ARGS[2])
rep = parse(Int, ARGS[3])
nx = parse(Int, ARGS[4])
ny = parse(Int, ARGS[5])

using KLIEPInference
using JLD
using LinearAlgebra
using SparseArrays

using ProximalBase, CoordinateDescent
include("/home/bkim6/code/refit_to_supp.jl")

filepath = "/home/bkim6/ising/"
file = jldopen(string(filepath, "$(graphtype)/$(graphtype)_rep_$(rep).jld"), "r")

X = read(file, "X")
Y = read(file, "Y")
close(file)

X = X[:, 1:nx]
Y = Y[:, 1:ny]

Ψx = Ψising(X .== 1.0)
Ψy = Ψising(Y .== 1.0)

p = size(Ψx, 1)

λ1 = round.(exp.(range(log(sqrt(2. * log(p) / nx)), stop=log(2. * sqrt(log(p) / nx)), length=5)), digits=3)
λ2 = sqrt(2. * log(p) / ny)

# rows : one-step vs double selection
# cols : λ1
θhat = zeros(2, 5)
σhat = zeros(2, 5)

for t = 1:5
    @show λ1[t]

    θtemp = spKLIEP(Ψx, Ψy, λ1[t], CD_KLIEP())
    θtemp = convert(SparseVector, θtemp)
    KLIEP_refit!(θtemp, Ψx, Ψy, union(idx, θtemp.nzind))

    H = KLIEP_Hessian(θtemp, Ψy)
    ω = Hinv_row(H, idx, λ2)
    ω = convert(SparseVector, ω)
    Hinv_row_refit!(ω, H, idx, ω.nzind)

    θhat[1, t] = KLIEP_debias(idx, θtemp, ω, Ψx, Ψy)
    σhat[1, t] = KLIEP_var(Ψx, Ψy, θtemp, ω)

    supp = union(idx, θtemp.nzind, ω.nzind)
    KLIEP_refit!(θtemp, Ψx, Ψy, supp)
    θhat[2, t] = θtemp[idx]
    δ = supp .== idx
    H = KLIEP_Hessian(θtemp[supp], Ψy[supp,:])
    ω = H\δ
    σhat[2, t] = KLIEP_var(Ψx[supp,:], Ψy[supp,:], θtemp[supp], ω)
end

@show θhat

@save "/home/bkim6/single/$(graphtype)/exp3_y_rep_$(rep).jld" θhat σhat
