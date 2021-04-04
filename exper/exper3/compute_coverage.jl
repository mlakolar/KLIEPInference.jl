using KLIEPInference
using JLD

scratch_dir = ARGS[1]

p = vcat(0.05:0.05:0.95, 0.90:0.01:0.99)

for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        coverage = zeros(length(p), 6)

        nrep = 0
        for rep in 1:1000
            try
                file = jldopen("$(scratch_dir)/res_$(m)_$(sgn)_$(rep).jld", "r")
                boot1 = read(file, "boot1")
                boot2 = read(file, "boot2")
                close(file)

                q = boot_quantile(boot1, p)
                covered1 = maximum(abs.(boot1.θhat)) .<= q

                q = boot_quantile(boot2, p)
                covered2 = maximum(abs.(boot2.θhat)) .<= q

                coverage[:, 1] .+= covered1
                coverage[:, 2] .+= covered2
                coverage[:, 3] .+= covered1 .== covered2

                q, w = boot_quantile_studentized(boot1, p)
                covered1 = maximum(abs.(boot1.θhat ./ w)) .<= q

                q, w = boot_quantile_studentized(boot2, p)
                covered2 = maximum(abs.(boot2.θhat ./ w)) .<= q

                coverage[:, 4] .+= covered1
                coverage[:, 5] .+= covered2
                coverage[:, 6] .+= covered1 .== covered2

                nrep += 1
            catch
                @warn "could not open $(m) $(sgn) $(rep)"
            end
        end

        coverage = broadcast(/, coverage, nrep)

        @save "coverage_$(m)_$(sgn).jld" coverage
    end
end
