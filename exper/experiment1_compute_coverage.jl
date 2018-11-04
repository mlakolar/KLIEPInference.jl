using KLIEPInference
using JLD

NUM_REP = 1000
coverage = 0
coverage_stud = 0

p = 10
sgn = 0
nx = 500
ny = 500

for rep=1:NUM_REP

    fname = "/scratch/midway2/mkolar/KLIEP/exp1/res_$(p)_$(sgn)_$(nx)_$(ny)_$(rep).jld"

    file = jldopen(fname, "r")
    res = read(file, "res")
    close(file)

    CI = simulCI(res)
    coverage += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0

    CI = simulCI(simulCIstudentized)
    coverage_stud += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 1 : 0
end

@show coverage / NUM_REP
@show coverage_stud / NUM_REP
