using KLIEPInference
using JLD

power = 0
power_stud = 0
not_done = 0

p = parse(Int,ARGS[1])
sgn = parse(Int, ARGS[2])
numChanges = parse(Int,ARGS[3])
lbInd = parse(Int,ARGS[4])
nx = parse(Int,ARGS[5])
ny = parse(Int,ARGS[6])
resPath = ARGS[7]
alpha = parse(Float64, ARGS[8])

NUM_REP = 1000

for rep=1:NUM_REP
    global power
    global not_done

    fname = "$(resPath)/res_$(p)_$(sgn)_$(numChanges)_$(lbInd)_$(nx)_$(ny)_$(rep).jld"

    try
        file = jldopen(fname, "r")
        res = read(file, "res")
        close(file)

        CI = simulCI(res, alpha)
        power += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 0 : 1
    catch
        not_done -= 1
    end
end

@show power / (NUM_REP + not_done)

not_done = 0

for rep=1:NUM_REP
    global power_stud
    global not_done

    fname = "$(resPath)/res_$(p)_$(sgn)_$(numChanges)_$(lbInd)_$(nx)_$(ny)_$(rep).jld"

    try
        file = jldopen(fname, "r")
        res = read(file, "res")
        close(file)

        CI = simulCIstudentized(res, alpha)
        power_stud += all(0 .<= CI[:, 2]) * all(0 .>= CI[:, 1]) ? 0 : 1
    catch
        not_done -= 1
    end
end

@show power_stud / (NUM_REP + not_done)
