### Steps to reproduce experimental results:

1. generate graphs: `./generate_graphs.sh`
2. run differential graph inference for 1000 repetitions: `./batch_exp3.sh`, may take a while. Change `-j NUMBER_OF_CORES` parameter fed to `parallel` to speed it up. 

   To run it on a slurm server: `./start_sbatch.sh`. Do not forget to set up the scratch folder by specifying `SCRATCH_DIR`.
3. compute coverage: `julia compute_coverage.jl "YOUR_SCRATCH_FOLDER"`.
4. plot result: `julia plot_res.jl` and the plot is generated as `res/*.png`.
