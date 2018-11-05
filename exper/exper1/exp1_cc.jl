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
sgn = parse(Int, ARGS[2])
nx = parse(Int,ARGS[3])
ny = parse(Int,ARGS[4])
resPath = ARGS[5]
alpha = parse(Float64, ARGS[6])

# resPath = "/scratch/midway2/mkolar/KLIEP/exp1/oracle"
NUM_REP = 1000



for rep=1:NUM_REP
    global coverage
    global coverage_stud
    global not_done

    fname = "$(resPath)/res_$(p)_$(sgn)_$(nx)_$(ny)_$(rep).jld"

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
