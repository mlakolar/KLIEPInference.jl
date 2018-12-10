using KLIEPInference
using JLD



p = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
numChanges = parse(Int,ARGS[3])
lbInd = parse(Int,ARGS[4])
nx = parse(Int,ARGS[5])
ny = parse(Int,ARGS[6])
resPath = ARGS[7]
outPrefix = ARGS[8]

# resPath = "/scratch/midway2/mkolar/KLIEP/exp2/oracle"
NUM_REP = 1000

results = Vector{BootstrapEstimates}(undef, NUM_REP)
for rep=1:NUM_REP
    global results
    fname = "$(resPath)/res_$(p)_$(sgn)_$(numChanges)_$(lbInd)_$(nx)_$(ny)_$(rep).jld"

    try
      file = jldopen(fname, "r")
      res = read(file, "res")
      close(file)

      results[rep] = deepcopy(res)
    catch
        @show fname
    end
end

@save "$(outPrefix)_p_$(p)_sgn_$(sgn)_numChanges_$(numChanges)_lbInd_$(lbInd)_nx_$(nx)_ny_$(ny).jld" results
                    