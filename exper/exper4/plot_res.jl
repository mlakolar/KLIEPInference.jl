using Distributions, StatsBase
using JLD, PyPlot

for m in 25
    for sgn in 1
        for numChanges in [1, 3, 5]
            file = jldopen("res/power_$(m)_$(sgn)_$(numChanges).jld", "r")
            power = read(file, "power")
            close(file)

            fig = figure(figsize=(3, 3), dpi=300)

            ax = subplot(1, 1, 1)

            scatter(0.0:0.05:0.5, power[1, :], marker="o")
            scatter(0.0:0.05:0.5, power[2, :], marker="v")

            minorticks_on()
            grid()
            grid(which="minor", ls="dotted", lw=".5")

            xlabel("δ", size="small")
            ylabel("proportion of rejections", size="small")
            ylim(0, 1.05)

            ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

            tight_layout()

            savefig("power_$(m)_$(sgn)_$(numChanges).png")

            close(fig)
        end
    end
end

for m in 25
    for sgn in 1
        for numChanges in [1, 3, 5]
            file = jldopen("res/power_$(m)_$(sgn)_$(numChanges).jld", "r")
            power = read(file, "power")
            close(file)

            fig = figure(figsize=(3, 3), dpi=300)

            ax = subplot(1, 1, 1)

            scatter(0.0:0.05:0.5, power[3, :], marker="o")
            scatter(0.0:0.05:0.5, power[4, :], marker="v")

            minorticks_on()
            grid()
            grid(which="minor", ls="dotted", lw=".5")

            xlabel("δ", size="small")
            ylabel("proportion of rejections", size="small")
            ylim(0, 1.05)

            ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

            tight_layout()

            savefig("power_studentized_$(m)_$(sgn)_$(numChanges).png")

            close(fig)
        end
    end
end
