using Distributions, StatsBase
using JLD, PyPlot

function qqplot(z; color="black")
    n = length(z)

    grid = [(1 / (n + 1)):(1 / (n + 1)):(1.0 - (1 / (n + 1)));]

    qz = quantile(z, grid)
    qd = quantile.(Ref(Distributions.Normal()), grid)

    lims = 3.290
    x = range(-lims, stop=lims)

    plot(x, x, color="grey", linestyle=":", linewidth=.25)
    scatter(qz, qd, s=.75, color=color)

    xlim([-lims, lims])
    ylim([-lims, lims])

    nothing
end

for m in [25, 50]
    for graph in ["chain1", "chain2", "tree1", "tree2"]
        file = jldopen("res/res_$(graph)_$(m).jld", "r")
        res = read(file, "res")

        x = range(-3.290, stop=3.290, length=100)
        y = pdf.(Normal(), x)

        fig = figure(figsize=(6,4), dpi=300)

        ax = subplot(2, 3, 1)
        qqplot(res[:, 1, 2], color="grey")
        qqplot(res[:, 2, 2])
        title("naive re-fitting", size="small")
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

        ax = subplot(2, 3, 2)
        qqplot(res[:, 1, 2], color="grey")
        qqplot(res[:, 3, 2])
        title("SparKLIE+1", size="small")
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

        ax = subplot(2, 3, 3)
        qqplot(res[:, 1, 2], color="grey")
        qqplot(res[:, 4, 2])
        title("SparKLIE+2", size="small")
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

        ax = subplot(2, 3, 4)
        plt[:hist](res[:, 2, 2], 100, density=true)
        plot(x,y)
        xlim(-3.290,3.290)
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

        ax = subplot(2, 3, 5)
        plt[:hist](res[:, 3, 2], 100, density=true)
        plot(x,y)
        xlim(-3.290,3.290)
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

        ax = subplot(2, 3, 6)
        plt[:hist](res[:, 4, 2], 100, density=true)
        plot(x,y)
        xlim(-3.290,3.290)
        ax[:tick_params]("both", labelsize="xx-small", length=2, pad=2)

        tight_layout()

        savefig("exper1_$(graph)_$(m).png")
        close(fig)
    end
end
