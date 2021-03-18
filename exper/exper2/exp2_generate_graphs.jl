using KLIEPInference
using JLD

file = jldopen("../exper1/graphs/chain1_25.jld", "r")
γy = pack(read(file, "Θy"))
close(file)

# setting 1. (none)
idx0 = KLIEPInference.trimap(5, 6)
for δ in -0.75:0.05:0.75
    γx = γy
    γx[idx0] += δ
    @save "graphs/params_exp2_set1_$(δ).jld" γx γy
end

# setting 2. (strong)
idx1 = KLIEPInference.trimap(4, 5)
idx2 = KLIEPInference.trimap(6, 7)
for δ in -0.75:0.05:0.75
    γx = γy
    γx[idx0] += δ
    γx[idx1] += 0.4
    γx[idx2] -= 0.4
    @save "graphs/params_exp2_set2_$(δ).jld" γx γy
end

# setting 3. (weak)
idx3 = KLIEPInference.trimap(4, 6)
idx4 = KLIEPInference.trimap(5, 7)
for δ in -0.75:0.05:0.75
    γx = γy
    γx[idx0] += δ
    γx[idx3] += 0.2
    γx[idx4] += 0.2
    @save "graphs/params_exp2_set3_$(δ).jld" γx γy
end

# setting 4. (mixed)
for δ in -0.75:0.05:0.75
    γx = γy
    γx[idx0] += δ
    γx[idx1] += 0.4
    γx[idx2] -= 0.4
    γx[idx3] += 0.2
    γx[idx4] += 0.2
    @save "graphs/params_exp2_set4_$(δ).jld" γx γy
end
