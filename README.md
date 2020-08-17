# glmGamPoi-Paper

> glmGamPoi: Fitting Gamma-Poisson Generalized Linear Models on Single Cell Count Data  
> Constantin Ahlmann-Eltze, Wolfgang Huber  
> bioRxiv 2020.08.13.249623; doi: https://doi.org/10.1101/2020.08.13.249623

This repository contains the code that was used to generate the figures for the paper that describes the [glmGamPoi](https://github.com/const-ae/glmGamPoi) package.

The scripts that call `glmGamPoi`, `DESeq2`, and `edgeR` are in the [benchmarks](https://github.com/const-ae/glmGamPoi-Paper/tree/master/benchmarks) folder:

* [run_benchmark-brain20k.R](https://github.com/const-ae/glmGamPoi-Paper/blob/master/benchmarks/run_benchmark-brain20k.R)
* [run_benchmark-mousegastrulation.R](https://github.com/const-ae/glmGamPoi-Paper/blob/master/benchmarks/run_benchmark-mousegastrulation.R)
* [run_benchmark-pbmc4k.R](https://github.com/const-ae/glmGamPoi-Paper/blob/master/benchmarks/run_benchmark-pbmc4k.R)
* [run_benchmark-pbmc68k.R](https://github.com/const-ae/glmGamPoi-Paper/blob/master/benchmarks/run_benchmark-pbmc68k.R)
* [run_asymptotics_benchmark.R](https://github.com/const-ae/glmGamPoi-Paper/blob/master/benchmarks/run_asymptotics_benchmark.R)
* [run_asymptotics_parameters_benchmark.R](https://github.com/const-ae/glmGamPoi-Paper/blob/master/benchmarks/run_asymptotics_parameters_benchmark.R)

Each script was run 5 times on a cluster and the results were extracted from the log files with the `gather_XXX_results.Rmd` scripts which are in the `benchmarks` folder as well.

The [data](https://github.com/const-ae/glmGamPoi-Paper/tree/master/data) folder contains the cleaned up tables that were used to produce the plots in the paper. 

The plots were generated with the [make_plots.Rmd](http://htmlpreview.github.io/?https://github.com/const-ae/glmGamPoi-Paper/blob/master/make_plots.nb.html) script.



