using JLD
using Printf

scratch_dir = ARGS[1]

for set in ["set1", "set2", "set3", "set4"]
    power = zeros(2, length(-0.75:0.05:0.75))
    
    i = 0
    for del in -0.75:0.05:0.75
        i += 1
        nrep = 0
        for batch in 1:10
            try
                file = jldopen("$(scratch_dir)/res_$(set)_$(del)_$(batch).jld", "r")
                T1 = read(file, "T1")
                T2 = read(file, "T2")
                close(file)

                power[1, i] += sum(T1)
                power[2, i] += sum(T2)

                nrep += 1
            catch
                @warn "could not open $(set) $(del) $(batch)"
            end
        end
	    power[1, i] = power[1, i] / (nrep * 100)
	    power[2, i] = power[2, i] / (nrep * 100)
    end

    @save "power_$(set).jld" power
end
