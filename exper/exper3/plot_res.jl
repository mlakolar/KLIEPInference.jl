using JLD, PyPlot

for m in [25, 50]
    fig = figure(figsize=(6,2), dpi=300)

    file = jldopen("res/coverage_$(m)_-1.jld", "r")
    coverage = read(file, "coverage")
    close(file)

    ax = subplot(1, 3, 1)
    plot(0.95:-0.05:0.05, 0.95:-0.05:0.05, color="grey", linestyle=":", linewidth=.25)
    plot(0.95:-0.05:0.05, coverage[1, :], marker="o")
    plot(0.95:-0.05:0.05, coverage[2, :], marker="o")
    title("sgn = -1", size="small")
    ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

    file = jldopen("res/coverage_$(m)_0.jld", "r")
    coverage = read(file, "coverage")
    close(file)

    ax = subplot(1, 3, 2)
    plot(0.95:-0.05:0.05, 0.95:-0.05:0.05, color="grey", linestyle=":", linewidth=.25)
    plot(0.95:-0.05:0.05, coverage[1, :], marker="o")
    plot(0.95:-0.05:0.05, coverage[2, :], marker="o")
    title("sgn = 0", size="small")
    ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

    file = jldopen("res/coverage_$(m)_1.jld", "r")
    coverage = read(file, "coverage")
    close(file)

    ax = subplot(1, 3, 3)
    plot(0.95:-0.05:0.05, 0.95:-0.05:0.05, color="grey", linestyle=":", linewidth=.25)
    plot(0.95:-0.05:0.05, coverage[1, :], marker="o")
    plot(0.95:-0.05:0.05, coverage[2, :], marker="o")
    title("sgn = +1", size="small")
    ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

    tight_layout()

    savefig("exper3_coverage_$(m).png")
    close(fig)

    fig = figure(figsize=(6,2), dpi=300)

    file = jldopen("res/coverage_$(m)_-1.jld", "r")
    coverage = read(file, "coverage")
    close(file)

    ax = subplot(1, 3, 1)
    plot(coverage[3, :], marker="o")
    title("sgn = -1", size="small")
    ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

    file = jldopen("res/coverage_$(m)_0.jld", "r")
    coverage = read(file, "coverage")
    close(file)

    ax = subplot(1, 3, 2)
    plot(coverage[3, :], marker="o")
    title("sgn = 0", size="small")
    ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

    file = jldopen("res/coverage_$(m)_1.jld", "r")
    coverage = read(file, "coverage")
    close(file)

    ax = subplot(1, 3, 3)
    plot(coverage[3, :], marker="o")
    title("sgn = +1", size="small")
    ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

    tight_layout()

    savefig("exper3_comparison_$(m).png")
    close(fig)
end
