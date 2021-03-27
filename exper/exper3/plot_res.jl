using JLD, PyPlot

for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        file = jldopen("res/coverage_$(m)_$(sgn).jld", "r")
        coverage = read(file, "coverage")
        close(file)

        # from 5% to 95%
        fig = figure(figsize=(3, 3), dpi=300)

        plot(0.0:0.05:1.0, 0.0:0.05:1.0, color="grey", linestyle=":", linewidth=.25)

        plot(0.05:0.05:0.95, coverage[1:19, 1], linewidth=.25, marker="o", markersize=2)
        plot(0.05:0.05:0.95, coverage[1:19, 2], linewidth=.25, marker="v", markersize=2)

        xlabel("target coverage", size="xx-small")
        ylabel("actual coverage", size="xx-small")

        xlim(0.0, 1.0)
        ylim(0.0, 1.0)

        tick_params("both", labelsize="xx-small")

        tight_layout()

        savefig("res/exper3_coverage_all_$(m)_$(sgn).png")

        close(fig)

        # from 5% to 95%, studentized
        fig = figure(figsize=(3, 3), dpi=300)

        plot(0.0:0.05:1.0, 0.0:0.05:1.0, color="grey", linestyle=":", linewidth=.25)

        plot(0.05:0.05:0.95, coverage[1:19, 4], linewidth=.25, marker="o", markersize=2)
        plot(0.05:0.05:0.95, coverage[1:19, 5], linewidth=.25, marker="v", markersize=2)

        xlabel("target coverage", size="xx-small")
        ylabel("actual coverage", size="xx-small")

        xlim(0.0, 1.0)
        ylim(0.0, 1.0)

        tick_params("both", labelsize="xx-small")

        tight_layout()

        savefig("res/exper3_coverage_all_studentized_$(m)_$(sgn).png")

        close(fig)

        # compare
        fig = figure(figsize=(3, 3), dpi=300)

        plot(0.05:0.05:0.95, coverage[1:19, 3], linewidth=.25, marker="o", markersize=2)

        xlabel("target coverage", size="xx-small")
        ylabel("proportion SparKLIE+1 == SparKLIE+2", size="xx-small")

        xlim(0.00, 1.0)
        ylim(0.88, 1.0)

        tick_params("both", labelsize="xx-small")

        tight_layout()

        savefig("res/exper3_compare_all_$(m)_$(sgn).png")

        close(fig)

        # compare, studentized
        fig = figure(figsize=(3, 3), dpi=300)

        plot(0.05:0.05:0.95, coverage[1:19, 6], linewidth=.25, marker="o", markersize=2)

        xlabel("target coverage", size="xx-small")
        ylabel("proportion SparKLIE+1 == SparKLIE+2", size="xx-small")

        xlim(0.00, 1.0)
        ylim(0.88, 1.0)

        tick_params("both", labelsize="xx-small")

        tight_layout()

        savefig("res/exper3_compare_all_studentized_$(m)_$(sgn).png")

        close(fig)

        # from 90% to 99%
        fig = figure(figsize=(3, 3), dpi=300)

        plot(0.89:0.01:1.00, 0.89:0.01:1.00, color="grey", linestyle=":", linewidth=.25)

        plot(0.90:0.01:0.99, coverage[20:end, 1], linewidth=.25, marker="o", markersize=2)
        plot(0.90:0.01:0.99, coverage[20:end, 2], linewidth=.25, marker="v", markersize=2)

        xlabel("target coverage", size="xx-small")
        ylabel("actual coverage", size="xx-small")

        xlim(0.89, 1.0)
        ylim(0.79, 1.0)

        tick_params("both", labelsize="xx-small")

        tight_layout()

        savefig("res/exper3_coverage_$(m)_$(sgn).png")

        close(fig)

        # from 90% to 99%, studentized
        fig = figure(figsize=(3, 3), dpi=300)

        plot(0.89:0.01:1.00, 0.89:0.01:1.00, color="grey", linestyle=":", linewidth=.25)

        plot(0.90:0.01:0.99, coverage[20:end, 4], linewidth=.25, marker="o", markersize=2)
        plot(0.90:0.01:0.99, coverage[20:end, 5], linewidth=.25, marker="v", markersize=2)

        xlabel("target coverage", size="xx-small")
        ylabel("actual coverage", size="xx-small")

        xlim(0.89, 1.0)
        ylim(0.79, 1.0)

        tick_params("both", labelsize="xx-small")

        tight_layout()

        savefig("res/exper3_coverage_studentized_$(m)_$(sgn).png")

        close(fig)
    end
end
