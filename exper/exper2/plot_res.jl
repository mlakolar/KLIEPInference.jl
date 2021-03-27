using Distributions, StatsBase
using JLD, PyPlot

for set in ["set1", "set2", "set3", "set4"]
    file = jldopen("res/power_$(set).jld", "r")
    power = read(file, "power")
    close(file)

    fig = figure(figsize=(3, 3), dpi=300)

    scatter(-.75:.05:.75, power[1, :], s=5, marker="o")
    scatter(-.75:.05:.75, power[2, :], s=5, marker="v")

    minorticks_on()
    grid()
    grid(which="minor", ls="dotted", lw=".5")

    xlabel("δ", size="xx-small")
    ylabel("proportion of rejections", size="xx-small")

    ylim(0, 1.05)

    tick_params("both", labelsize="xx-small")

    tight_layout()

    savefig("power_$(set).png")

    close(fig)
end
