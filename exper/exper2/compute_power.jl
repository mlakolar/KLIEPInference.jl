using JLD

for set in ["set1", "set2", "set3", "set4"]
    power = zeros(2, length(-0.75:0.05:0.75))
    i = 0
    for del in -0.75:0.05:0.75
        i += 1
        nrep = 0
        for rep in 1:1000
            try
                file = jldopen("/scratch/midway2/byolkim/exper2/res_$(set)_$(del)_$(rep).jld", "r")
                T1 = read(file, "T1")
                T2 = read(file, "T2")
                close(file)

                power[1, i] += T1
                power[2, i] += T2

                nrep += 1
            catch
                @warn "Could not open rep $(rep)"
            end
        end
    end
    power = broadcast(/, power, nrep)

    @save "/scratch/midway2/byolkim/exper2/power_$(set).jld" power
end
