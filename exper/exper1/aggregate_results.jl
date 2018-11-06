using KLIEPInference
using JLD

p = parse(Int,ARGS[1])
sgn = parse(Int, ARGS[2])
nx = parse(Int,ARGS[3])
ny = parse(Int,ARGS[4])
resPath = ARGS[5]
outPrefix = ARGS[6]

# resPath = "/scratch/midway2/mkolar/KLIEP/exp1/oracle"
NUM_REP = 1000

results = Vector{BootstrapEstimates}(undef, NUM_REP)
for rep=1:NUM_REP
    global results
    fname = "$(resPath)/res_$(p)_$(sgn)_$(nx)_$(ny)_$(rep).jld"

    try
      file = jldopen(fname, "r")
      res = read(file, "res")
      close(file)

      results[rep] = deepcopy(res)
    catch
        @show fname
    end
end

@save "$(outPrefix)_p_$(p)_sgn_$(sgn)_nx_$(nx)_ny_$(ny).jld" results
