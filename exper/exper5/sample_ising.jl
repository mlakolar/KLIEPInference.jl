graphtype = ARGS[1]
rep = parse(Int,ARGS[2])
nx = parse(Int,ARGS[3])
ny = parse(Int,ARGS[4])

if isfile("/home/bkim6/balanced/$(graphtype)/$(graphtype)_rep_$(rep).jld")
    exit()
end

using KLIEPInference
using Distributions
using Random
using SparseArrays
using JLD

import KLIEPInference.trimap

file = jldopen("/home/bkim6/graphs/$(graphtype).jld", "r")
Θx = read(file, "Θx")
Θy = read(file, "Θy")
close(file)

Θx = sparse(Θx)
Θy = sparse(Θy)

m = size(Θx, 1)
p = div(m * (m-1), 2)
θx = zeros(p)
θy = zeros(p)
for (i,j) in zip(findnz(Θx)...)
    θx[trimap(i,j)] = Θx[i,j]
    Θy[trimap(i,j)] = Θy[i,j]
end

Random.seed!(123 + rep)
spl = IsingSampler(θx; thin=2000)
X = rand(spl, nx)
spl = IsingSampler(θy; thin=2000)
Y = rand(spl, ny)

if ~ispath("/home/bkim6/balanced/$(graphtype)/")
    mkpath("/home/bkim6/balanced/$(graphtype)/")
end
@save "/home/bkim6/balanced/$(graphtype)/$(graphtype)_rep_$(rep).jld" X Y
