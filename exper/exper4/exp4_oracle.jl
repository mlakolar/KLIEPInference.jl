# for the balanced sample sizes

graphdir = ARGS[1]
graphtype = ARGS[2]
idx = parse(Int, ARGS[3])
rep = parse(Int, ARGS[4])
nx = parse(Int, ARGS[5])
ny = parse(Int, ARGS[6])

using KLIEPInference
using JLD
using LinearAlgebra
using SparseArrays

file = jldopen("/home/bkim6/graphs/$(graphtype).jld", "r")
Θx = read(file, "Θx")
Θy = read(file, "Θy")
close(file)

Θ = Θx - Θy
θ = Θ[findall(triu(ones(Bool,size(Θ)),1))]
θ = sparse(θ)

supp = θ.nzind
δ = supp .== idx

file = jldopen("$(graphdir)/$(graphtype)/$(graphtype)_rep_$(rep).jld", "r")
X = read(file, "X")
Y = read(file, "Y")
close(file)

X = X[:, 1:nx]
Y = Y[:, 1:ny]

Ψx = Ψising(X .== 1.)
Ψy = Ψising(Y .== 1.)

Ψx = Ψx[supp, :]
Ψy = Ψy[supp, :]

θhat = zeros(3)
σhat = zeros(3)

# Base y
θtemp = KLIEP(Ψx, Ψy, CD_KLIEP())
θhat[1] = θtemp[δ][1]

H = KLIEP_Hessian(θtemp, Ψy)
ω = H\δ
σhat[1] = KLIEP_var(Ψx, Ψy, θtemp, ω)

# Base x
θtemp = KLIEP(Ψy, Ψx, CD_KLIEP())
θhat[2] = -θtemp[δ][1]

H = KLIEP_Hessian(θtemp, Ψx)
ω = H\δ
σhat[2] = KLIEP_var(Ψy, Ψx, θtemp, ω)

# SymKLIEP
θtemp = SymKLIEP(Ψx, Ψy, CD_SymKLIEP())
θhat[3] = θtemp[δ][1]

H = SymKLIEP_Hessian(θtemp, Ψx, Ψy)
ω = H\δ
σhat[3] = SymKLIEP_var(Ψx, Ψy, θtemp, ω)

@show θhat

if ~ispath("/home/bkim6/oracle_balanced/$(graphtype)/")
    mkpath("/home/bkim6/oracle_balanced/$(graphtype)/")
end
@save "/home/bkim6/oracle_balanced/$(graphtype)/$(graphtype)_rep_$(rep).jld" θhat σhat
