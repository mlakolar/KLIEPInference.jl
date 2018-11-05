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
NUM_REP = parse(Int,ARGS[5])

file = jldopen("params_exp2_$(p)_$(sgn).jld", "r")
θx = read(file, "θx")
θy = read(file, "θy")
close(file)

Δ = θx - θy

for rep=1:NUM_REP
    global coverage
    global coverage_stud
    global not_done

    fname = "/scratch/midway2/mkolar/KLIEP/exp2/res_$(p)_$(sgn)_$(nx)_$(ny)_$(rep).jld"

    try 
      file = jldopen(fname, "r")
      res = read(file, "res")
      close(file)

      CI = simulCI(res)
      # coverage += all(Δ .<= CI[:, 2]) * all(Δ .>= CI[:, 1]) ? 1 : 0
      coverage += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

      CI = simulCIstudentized(res)
      # coverage_stud += all(Δ .<= CI[:, 2]) * all(Δ .>= CI[:, 1]) ? 1 : 0
      coverage_stud += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0
    catch
      not_done -= 1
    end
end

@show coverage / (NUM_REP + not_done)
@show coverage_stud / (NUM_REP + not_done)

