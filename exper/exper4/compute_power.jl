using KLIEPInference
using JLD

scratch_dir = ARGS[1]

p = 0.05:0.05:0.95

for m in 25
    for sgn in 1
        for numChanges in [1, 3, 5]
            println("computing power for $(m) $(sgn) $(numChanges) ...")

            power = zeros(11, 4)
            for lbInd = 1:11
                file = jldopen("graphs/params_exp4_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld", "r")
                γx = read(file, "γx")
                γy = read(file, "γy")
                θ = γx - γy
                close(file)

                coverage = zeros(length(p), 4)

                nrep = 0
                for rep in 1:1000
                    try
                        file = jldopen("$(scratch_dir)/res_$(m)_$(sgn)_$(numChanges)_$(lbInd)_$(rep).jld", "r")
                        boot1 = read(file, "boot1")
                        boot2 = read(file, "boot2")
                        close(file)

                        # coverage
                        q = boot_quantile(boot1, p)
                        coverage[:, 1] .+= maximum(abs.(boot1.θhat - θ)) .<= q

                        q = boot_quantile(boot2, p)
                        coverage[:, 2] .+=  maximum(abs.(boot2.θhat - θ)) .<= q

                        q, w = boot_quantile_studentized(boot1, p)
                        coverage[:, 3] .+= maximum(abs.((boot1.θhat - θ) ./ w)) .<= q

                        q, w = boot_quantile_studentized(boot2, p)
                        coverage[:, 4] .+= maximum(abs.((boot2.θhat - θ) ./ w)) .<= q

                        # power
                        q = boot_quantile(boot1, 0.95)
                        power[lbInd, 1] += maximum(abs.(boot1.θhat)) > q

                        q = boot_quantile(boot2, 0.95)
                        power[lbInd, 2] += maximum(abs.(boot2.θhat)) > q

                        q, w = boot_quantile_studentized(boot1, 0.95)
                        power[lbInd, 3] += maximum(abs.(boot1.θhat ./ w)) > q

                        q, w = boot_quantile_studentized(boot2, 0.95)
                        power[lbInd, 4] += maximum(abs.(boot2.θhat ./ w)) > q

                        nrep += 1
                    catch
                        @warn "could not open $(m) $(sgn) $(numChanges) $(lbInd) $(rep)"
                    end
                end

                coverage = broadcast(/, coverage, nrep)

                @save "res/coverage_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld" coverage

                power[lbInd, :] ./= nrep
            end

            @save "res/power_$(m)_$(sgn)_$(numChanges).jld" power

            println("... done!")
        end
    end
end
