using KLIEPInference
using JLD

resPath = ARGS[1]
outPrefix = ARGS[2]

# resPath = "/scratch/midway2/mkolar/KLIEP/exp1/oracle"
NUM_REP = 1000

results = Vector{BootstrapEstimates}(undef, NUM_REP)
for rep=1:NUM_REP
    global results
    fname = "$(resPath)/res_$(rep).jld"

    try
      file = jldopen(fname, "r")
      res = read(file, "res")
      close(file)

      results[rep] = deepcopy(res)
    catch
        @show fname
    end
end

@save "$(outPrefix).jld" results
