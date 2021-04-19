using JLD, PyPlot

# Figure 1
fig = figure(figsize=(6, 6), dpi=300)

k = 0
for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        file = jldopen("res/coverage_$(m)_$(sgn).jld", "r")
        coverage = read(file, "coverage")
        close(file)

        global k += 1

        ax = subplot(3, 3, k)

        plot(0.0:0.05:1.0, 0.0:0.05:1.0, color="grey", linestyle=":", linewidth=.25)

        plot(0.05:0.05:0.95, coverage[1:19, 1], linewidth=.25, marker="o", markersize=2)
        plot(0.05:0.05:0.95, coverage[1:19, 2], linewidth=.25, marker="v", markersize=2)

        xlim(0.0, 1.0)
        ylim(0.0, 1.0)

        tick_params("both", labelsize="xx-small")

        title("m = $(m), sgn = $(sgn)", size="x-small")

        if m == 100 && sgn == 0
            xlabel(L"$1 - \alpha$", size="x-small")
        end
        if m == 50 && sgn == -1
            ylabel(L"proportion $T \leq \hat c_{T, 1-\alpha}$", size="x-small")
        end
    end
end

tight_layout()

savefig("res/exper3_coverage.png")

close(fig)

# studentized
fig = figure(figsize=(6, 6), dpi=300)

k = 0
for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        file = jldopen("res/coverage_$(m)_$(sgn).jld", "r")
        coverage = read(file, "coverage")
        close(file)

        global k += 1

        ax = subplot(3, 3, k)

        plot(0.0:0.05:1.0, 0.0:0.05:1.0, color="grey", linestyle=":", linewidth=.25)

        plot(0.05:0.05:0.95, coverage[1:19, 4], linewidth=.25, marker="o", markersize=2)
        plot(0.05:0.05:0.95, coverage[1:19, 5], linewidth=.25, marker="v", markersize=2)

        xlim(0.0, 1.0)
        ylim(0.0, 1.0)

        tick_params("both", labelsize="xx-small")

        title("m = $(m), sgn = $(sgn)", size="x-small")

        if m == 100 && sgn == 0
            xlabel(L"$1 - \alpha$", size="x-small")
        end
        if m == 50 && sgn == -1
            ylabel(L"proportion $T \leq \hat c_{T, 1-\alpha}$", size="x-small")
        end
    end
end

tight_layout()

savefig("res/exper3_coverage_studentized.png")

close(fig)

# compare SparKLIE+1 vs SparKLIE+2
fig = figure(figsize=(6, 6), dpi=300)

k = 0
for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        file = jldopen("res/coverage_$(m)_$(sgn).jld", "r")
        coverage = read(file, "coverage")
        close(file)

        global k += 1

        ax = subplot(3, 3, k)

        plot(0.05:0.05:0.95, coverage[1:19, 3], linewidth=.25, marker="o", markersize=2)

        xlim(0.00, 1.0)
        ylim(0.75, 1.0)

        tick_params("both", labelsize="xx-small")

        title("m = $(m), sgn = $(sgn)", size="x-small")

        if m == 100 && sgn == 0
            xlabel(L"$1 - \alpha$", size="x-small")
        end
        if m == 50 && sgn == -1
            ylabel("proportion SparKLIE+1 == SparKLIE+2", size="x-small")
        end
    end
end

tight_layout()

savefig("res/exper3_compare.png")

close(fig)

# compare SparKLIE+1 vs SparKLIE+2, studentized
fig = figure(figsize=(6, 6), dpi=300)

k = 0
for m in [25, 50, 100]
    for sgn in [-1, 0, 1]
        file = jldopen("res/coverage_$(m)_$(sgn).jld", "r")
        coverage = read(file, "coverage")
        close(file)

        global k += 1

        ax = subplot(3, 3, k)

        plot(0.05:0.05:0.95, coverage[1:19, 6], linewidth=.25, marker="o", markersize=2)

        xlim(0.00, 1.0)
        ylim(0.75, 1.0)

        tick_params("both", labelsize="xx-small")

        title("m = $(m), sgn = $(sgn)", size="x-small")

        if m == 100 && sgn == 0
            xlabel(L"$1 - \alpha$", size="x-small")
        end
        if m == 50 && sgn == -1
            ylabel("proportion SparKLIE+1 == SparKLIE+2", size="x-small")
        end
    end
end

tight_layout()

savefig("res/exper3_compare_studentized.png")

close(fig)
