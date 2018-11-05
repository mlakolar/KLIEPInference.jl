using KLIEPInference
using JLD

coverage = 0
coverage_stud = 0
not_done = 0

# p = 10
# sgn = 0
# nx = 500
# ny = 500


p = parse(Int,ARGS[1])
sgn = parse(Int,ARGS[2])
numChanges = parse(Int,ARGS[3])
lbInd = parse(Int,ARGS[4])
nx = parse(Int,ARGS[5])
ny = parse(Int,ARGS[6])
resPath = ARGS[7]
alpha = parse(Float64, ARGS[8])


# resPath = "/scratch/midway2/mkolar/KLIEP/exp2/oracle"
NUM_REP = 1000



for rep=1:NUM_REP
    global coverage
    global coverage_stud
    global not_done

    fname = "$(resPath)/res_$(p)_$(sgn)_$(numChanges)_$(lbInd)_$(nx)_$(ny)_$(rep).jld"
                        
    try 
      file = jldopen(fname, "r")
      res = read(file, "res")
      close(file)

      CI = simulCI(res, alpha)
      coverage += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

      CI = simulCIstudentized(res, alpha)
      coverage_stud += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0
    catch
      not_done -= 1
    end
end

@show coverage / (NUM_REP + not_done)
@show coverage_stud / (NUM_REP + not_done)
