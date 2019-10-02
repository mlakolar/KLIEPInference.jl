# Oracle comparisons

using KLIEPInference
using Distributions, StatsBase
using LinearAlgebra
using SparseArrays
using Random
using JLD

# using ProximalBase, CoordinateDescent
# include("/home/bkim6/code/utils_single_edge.jl")
# include("/home/bkim6/code/utils_SymKLIEPLoss.jl")
# include("/home/bkim6/code/solver_SymKLIEP.jl")

graphtype = ARGS[1]
rep = parse(Int,ARGS[2])
ind1 = parse(Int,ARGS[3])
nx = parse(Int,ARGS[4])
ny = parse(Int,ARGS[5])

file = jldopen("/home/bkim6/graphs/$(graphtype).jld", "r")
Θx = read(file, "Θx")
Θy = read(file, "Θy")
close(file)

Θ = Θx - Θy
θ = Θ[findall(triu(ones(Bool,size(Θ)),1))]
θ = sparse(θ)

supp = θ.nzind
ek = supp .== ind1

file = jldopen("/home/bkim6/ising/$(graphtype)/$(graphtype)_rep_$(rep).jld", "r")
X = read(file, "X")
Y = read(file, "Y")
close(file)

X = X[:, 1:nx]
Y = Y[:, 1:ny]

Ψx = Ψising(X .== 1.0)
Ψy = Ψising(Y .== 1.0)

Ψx = Ψx[supp, :]
Ψy = Ψy[supp, :]

println("keeping x and y the same as in the paper...")
θhat = KLIEP(Ψx, Ψy, CD_KLIEP())
θhat_y = θhat[ek][1]

H = KLIEP_Hessian(θhat, Ψy)
ω = H\ek
σhat_y = KLIEP_var(Ψx, Ψy, θhat, ω)
@show θhat_y σhat_y

println("switching x and y...")
θhat = KLIEP(Ψy, Ψx, CD_KLIEP())
θhat_x = -θhat[ek][1]

H = KLIEP_Hessian(θhat, Ψx)
ω = H\ek
σhat_x = KLIEP_var(Ψy, Ψx, θhat, ω)
@show θhat_x σhat_x

println("symmetrized procedure...")
θhat = SymKLIEP(Ψx, Ψy, CD_SymKLIEP())
θhat_sym = θhat[ek][1]

H = SymKLIEP_Hessian(θhat, Ψx, Ψy)
ω = H\ek
σhat_sym = SymKLIEP_var(Ψx, Ψy, θhat, ω)
@show θhat_sym σhat_sym

@save "/home/bkim6/oracle/$(graphtype)/$(graphtype)_rep_$(rep).jld" θhat_y σhat_y θhat_x σhat_x θhat_sym σhat_sym
