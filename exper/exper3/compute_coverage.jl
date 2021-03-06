using JLD

for m in [25, 50]
    for sgn in [-1, 0, 1]
        coverage = zeros(Int, 3, length(0.05:0.05:0.95));

        nrep = 0
        for rep in 1:1000
            try
                file = jldopen("/scratch/midway2/byolkim/exper3/res_$(m)_$(sgn)_$(rep).jld", "r")
                res = read(file, "res")
                close(file)

                coverage[1, :] .+= res[1, :]
                coverage[2, :] .+= res[2, :]
                coverage[3, :] .+= res[1, :] .=== res[2, :]

                nrep += 1
            catch
                @warn "Could not open rep $(rep)"
            end
        end

        coverage = broadcast(/, coverage, nrep)

        @save "/scratch/midway2/byolkim/exper3/coverage_$(m)_$(sgn).jld" coverage
    end
end
