using JLD, PyPlot

fig = figure(figsize=(6, 6), dpi=300)

k = 0
for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        global k += 1

        ax = subplot(3, 3, k)

        file = jldopen("res/coverage_$(m)_$(sgn).jld", "r")
        coverage = read(file, "coverage")
        close(file)

        plot(0.95:-0.05:0.05, 0.95:-0.05:0.05, color="grey", linestyle=":", linewidth=.25)

        plot(0.95:-0.05:0.05, coverage[1, :], linewidth=.25, marker="o", markersize=2)
        plot(0.95:-0.05:0.05, coverage[2, :], linewidth=.25, marker="v", markersize=2)

        if m === 25
            ax.set_title("sgn = $(sgn)", size="small")
        end

        if m === 100
            ax.set_xlabel("α", size="small")
        end

        if sgn === -1
            ax.set_ylabel("m = $(m)", size="small")
        end

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)
    end
end

tight_layout()

savefig("exper3_coverage.png")

close(fig)
