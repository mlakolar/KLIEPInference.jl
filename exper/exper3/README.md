### Steps to reproduce experimental results: 

1. generate graphs: `./generate_graphs.sh`
2. run differential graph inference for 1000 repetitions: `./batch_exp3.sh`, may take a while. 

   To run it on a slurm server: `./start_sbatch.sh`. 
3. compute coverage: `julia compute_coverage.jl`. 
4. plot result: `julia plot_res.jl` and the plot is generated as `exper3_coverage.png`.

![](exper3_coverage.png)