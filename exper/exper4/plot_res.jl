using DelimitedFiles, PyPlot

res, header = readdlm("res/exper4.csv", ',', Int64, '\n'; header=true)

for m in [25, 50, 100]
    fig = figure(figsize=(3.14961, 3.14961), dpi=1000)

    power = res[findall((res[:, 1] .== m) .& (res[:, 2] .== 1)), 4] ./ 1000

    scatter(0.0:0.05:0.5, power, s=5, marker="o")

    power = res[findall((res[:, 1] .== m) .& (res[:, 2] .== 3)), 4] ./ 1000

    scatter(0.0:0.05:0.5, power, s=5, marker="v")

    power = res[findall((res[:, 1] .== m) .& (res[:, 2] .== 5)), 4] ./ 1000

    scatter(0.0:0.05:0.5, power, s=5, marker="s")

    minorticks_on()
    grid()
    grid(which="minor", ls="dotted", lw=".5")

    title("without Studentization", size="xx-small")
    xlabel("δ", size="xx-small")
    ylabel("proportion of rejections", size="xx-small")

    ylim(0, 1.05)

    tick_params("both", labelsize="xx-small")

    tight_layout()

    savefig("res/exper4_$(m).pdf")

    close(fig)

    fig = figure(figsize=(3.14961, 3.14961), dpi=1000)

    power = res[findall((res[:, 1] .== m) .& (res[:, 2] .== 1)), 5] ./ 1000

    scatter(0.0:0.05:0.5, power, s=5, marker="o")

    power = res[findall((res[:, 1] .== m) .& (res[:, 2] .== 3)), 5] ./ 1000

    scatter(0.0:0.05:0.5, power, s=5, marker="v")

    power = res[findall((res[:, 1] .== m) .& (res[:, 2] .== 5)), 5] ./ 1000

    scatter(0.0:0.05:0.5, power, s=5, marker="s")

    minorticks_on()
    grid()
    grid(which="minor", ls="dotted", lw=".5")

    title("Studentized", size="xx-small")
    xlabel("δ", size="xx-small")
    ylabel("proportion of rejections", size="xx-small")

    ylim(0, 1.05)

    tick_params("both", labelsize="xx-small")

    tight_layout()

    savefig("res/exper4_$(m)_studentized.pdf")

    close(fig)
end
