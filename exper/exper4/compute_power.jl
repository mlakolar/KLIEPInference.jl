using JLD
using Printf

scratch_dir = ARGS[1]

for m in 25
    for sgn in 1
        for numChanges in [1, 3, 5]
            println("computing power for $(m) $(sgn) $(numChanges) ...")
            power = zeros(4, 11)
            for lbInd = 1:11
                nrep = 0
                for rep in 1:1000
                    try
                        file = jldopen("$(scratch_dir)/res_$(m)_$(sgn)_$(numChanges)_$(lbInd)_$(rep).jld", "r")
                        T1 = read(file, "T1")
                        T2 = read(file, "T2")
                        W1 = read(file, "W1")
                        W2 = read(file, "W2")
                        close(file)

                        power[1, lbInd] += T1
                        power[2, lbInd] += T2
                        power[3, lbInd] += W1
                        power[4, lbInd] += W2

                        nrep += 1
                    catch
                        @warn "Could not open rep $(m) $(sgn) $(numChanges) $(lbInd) $(rep)"
                    end
                end
                power[:, lbInd] ./= nrep
            end
            @save "power_$(m)_$(sgn)_$(numChanges).jld" power
            println("... done!")
        end
    end
end
