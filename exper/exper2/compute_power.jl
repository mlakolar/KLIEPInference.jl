using JLD
using Printf

scratch_dir = ARGS[1]

for set in ["set1", "set2", "set3", "set4"]
    power = zeros(2, length(-0.75:0.05:0.75))
    i = 0
    for del in -0.75:0.05:0.75
        i += 1
        nrep = 0
        for rep in 1:1000
            try
                file = jldopen("$(scratch_dir)/exp2_$(set)_$(@sprintf("%.2f", del))/res_$(set)_$(del)_$(rep).jld", "r")
                T1 = read(file, "T1")
                T2 = read(file, "T2")
                close(file)

                power[1, i] += T1
                power[2, i] += T2

                nrep += 1
            catch
                @warn "Could not open rep $(set) $(del) $(rep)"
            end
        end
	power[1, i] = power[1, i] / nrep
	power[2, i] = power[2, i] / nrep
    end


    @save "power_$(set).jld" power
end
