using JLD, PyPlot

for m in 25
    for sgn in 0
        for numChanges in [1, 3, 5]
            file = jldopen("res/power_$(m)_$(sgn)_$(numChanges).jld", "r")
            power = read(file, "power")
            close(file)

            fig = figure(figsize=(3, 3), dpi=300)

            scatter(0.0:0.05:0.4, power[1:9, 1], s=5, marker="o")
            scatter(0.0:0.05:0.4, power[1:9, 2], s=5, marker="v")

            minorticks_on()
            grid()
            grid(which="minor", ls="dotted", lw=".5")

            xlabel("δ", size="xx-small")
            ylabel("proportion of rejections", size="xx-small")

            ylim(0, 1.05)

            tick_params("both", labelsize="xx-small")

            tight_layout()

            savefig("res/power_$(m)_$(sgn)_$(numChanges).png")

            close(fig)
        end
    end
end

for m in 25
    for sgn in 0
        for numChanges in [1, 3, 5]
            file = jldopen("res/power_$(m)_$(sgn)_$(numChanges).jld", "r")
            power = read(file, "power")
            close(file)

            fig = figure(figsize=(3, 3), dpi=300)

            scatter(0.0:0.05:0.4, power[1:9, 3], s=5, marker="o")
            scatter(0.0:0.05:0.4, power[1:9, 4], s=5, marker="v")

            minorticks_on()
            grid()
            grid(which="minor", ls="dotted", lw=".5")

            xlabel("δ", size="xx-small")
            ylabel("proportion of rejections", size="xx-small")

            ylim(0, 1.05)

            tick_params("both", labelsize="xx-small")

            tight_layout()

            savefig("res/power_studentized_$(m)_$(sgn)_$(numChanges).png")

            close(fig)
        end
    end
end

for m in 25
    for sgn in 0
        fig = figure(figsize=(27, 9), dpi=300)

        k = 0
        for numChanges in [1, 3, 5]
            for lbInd = 1:9
                k += 1

                ax = subplot(3, 9, k)

                file = jldopen("res/coverage_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld", "r")
                coverage = read(file, "coverage")
                close(file)

                plot(0.0:0.05:1.0, 0.0:0.05:1.0, color="grey", linestyle=":", linewidth=.25)

                plot(0.05:0.05:0.95, coverage[1:19, 1], linewidth=.25, marker="o", markersize=2)
                plot(0.05:0.05:0.95, coverage[1:19, 2], linewidth=.25, marker="v", markersize=2)

                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0.0, 1.0)

                ax[:tick_params]("both", labelsize=2)
            end
        end

        tight_layout()

        savefig("res/coverage_$(m)_$(sgn).png")

        close(fig)
    end
end

for m in 25
    for sgn in 0
        fig = figure(figsize=(27, 9), dpi=300)

        k = 0
        for numChanges in [1, 3, 5]
            for lbInd = 1:9
                k += 1

                ax = subplot(3, 9, k)

                file = jldopen("res/coverage_$(m)_$(sgn)_$(numChanges)_$(lbInd).jld", "r")
                coverage = read(file, "coverage")
                close(file)

                plot(0.0:0.05:1.0, 0.0:0.05:1.0, color="grey", linestyle=":", linewidth=.25)

                plot(0.05:0.05:0.95, coverage[1:19, 3], linewidth=.25, marker="o", markersize=2)
                plot(0.05:0.05:0.95, coverage[1:19, 4], linewidth=.25, marker="v", markersize=2)

                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0.0, 1.0)

                ax[:tick_params]("both", labelsize=2)
            end
        end

        tight_layout()

        savefig("res/coverage_studentized_$(m)_$(sgn).png")

        close(fig)
    end
end
