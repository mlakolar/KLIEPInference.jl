# Repeat the single-edge experiment, but with x- and y- switched, and with symmetrized loss

using KLIEPInference
using Distributions
using LinearAlgebra
using SparseArrays
using Random
using JLD

graphtype = ARGS[1]
idx = parse(Int, ARGS[2])
rep = parse(Int, ARGS[3])
λ1 = parse(Float64, ARGS[4])
nx = parse(Int, ARGS[5])
ny = parse(Int, ARGS[6])


filepath = "/home/bkim6/ising/$(graphtype)/"
file = jldopen(string(filepath, "$(graphtype)_rep_$(rep).jld"), "r")

X = read(file, "X")
Y = read(file, "Y")
close(file)

X = X[:, 1:nx]
Y = Y[:, 1:ny]

Ψx = Ψising(X .== 1.0)
Ψy = Ψising(Y .== 1.0)

p = size(Ψx, 1)


# rows : one-step vs double selection
# cols : y-based, x-based, symmetric
θhat = zeros(2, 3)
σhat = zeros(2, 3)


println("nx < ny")
println("step 1")
@show λ1
θhaty = spKLIEP(Ψx, Ψy, λ1, CD_KLIEP())
spKLIEP_refit!(θhaty, Ψx, Ψy)
θhaty = SparseVector(θhaty)
println("step 1 done")

println("step 2")
@show λ2 = sqrt(2. * log(p) / ny)
H = KLIEP_Hessian(θhaty, Ψy)
ω = Hinv_row(H, idx, λ2)
ω = convert(SparseVector, ω)
println("step 2 done")

println("de-bias")
θhat[1, 1] = KLIEP_debias(idx, θhaty, ω, Ψx, Ψy)
σhat[1, 1] = KLIEP_var(Ψx, Ψy, θhaty, ω)

supp = union(θhaty.nzind, ω.nzind, idx)
ek = supp .== idx
θtilde = KLIEP(Ψx[supp, :], Ψy[supp, :], CD_KLIEP())
θhat[2, 1] = θtilde[ek][1]

H = KLIEP_Hessian(θhaty[supp], Ψy[supp,:])
ω = H\ek
σhat[2, 1] = KLIEP_var(Ψx[supp,:], Ψy[supp,:], θhaty[supp], ω)


println("nx > ny")
println("step 1")
@show λ1
θhatx = spKLIEP(Ψy, Ψx, λ1, CD_KLIEP())
spKLIEP_refit!(θhatx, Ψy, Ψx)
θhatx = SparseVector(θhatx)
println("step 1 done")

println("step 2")
@show λ2 = sqrt(2. * log(p) / nx)
H = KLIEP_Hessian(θhatx, Ψx)
ω = Hinv_row(H, idx, λ2)
ω = convert(SparseVector, ω)
println("step 2 done")

println("de-bias")
θhat[1, 2] = KLIEP_debias(idx, θhatx, ω, Ψy, Ψx)
σhat[1, 2] = KLIEP_var(Ψy, Ψx, θhatx, ω)

supp = union(θhatx.nzind, ω.nzind, idx)
ek = supp .== idx
θtilde = KLIEP(Ψy[supp, :], Ψx[supp, :], CD_KLIEP())
θhat[2, 2] = θtilde[ek][1]

H = KLIEP_Hessian(θhatx[supp], Ψx[supp,:])
ω = H\ek
σhat[2, 2] = KLIEP_var(Ψy[supp,:], Ψx[supp,:], θhatx[supp], ω)


println("symmetric")
println("step 1")
@show λ1 = 2 * λ1
θhatb = spSymKLIEP(Ψx, Ψy, λ1, CD_SymKLIEP())
spSymKLIEP_refit!(θhatb, Ψx, Ψy)
θhatb = SparseVector(θhatb)
println("step 1 done")

println("step 2")
@show λ2 = 0.5 * (sqrt(2. * log(p) / nx) + sqrt(2. * log(p) / ny))
H = SymKLIEP_Hessian(θhatb, Ψx, Ψy)
ω = Hinv_row(H, idx, λ2)
ω = convert(SparseVector, ω)
println("step 2 done")

println("de-bias")
θhat[1, 3] = SymKLIEP_debias(idx, θhatb, ω, Ψx, Ψy)
σhat[1, 3] = SymKLIEP_var(Ψx, Ψy, θhatb, ω)

supp = union(θhatb.nzind, ω.nzind, idx)
ek = supp .== idx
θtilde = SymKLIEP(Ψx[supp, :], Ψy[supp, :], CD_SymKLIEP())
θhat[2, 3] = θtilde[ek][1]

H = SymKLIEP_Hessian(θhatb[supp], Ψx[supp,:], Ψy[supp,:])
ω = H\ek
σhat[2, 3] = SymKLIEP_var(Ψx[supp,:], Ψy[supp,:], Array(θhatb[supp]), ω)

@save "/home/bkim6/single/$(graphtype)/$(graphtype)_rep_$(rep).jld" θhat σhat
