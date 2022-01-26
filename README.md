# BiTSC2
Bayesian inference of Tumor clonal Tree by joint analysis of Single-Cell SNV and CNA data

## Software dependencies

BiTSC2 is written with `R` and `C++`. Before implementing our software, please install the following packages in `R`:

data manipulation: `tidyr`, `reshape`, `dplyr`

Bayesian analysis: `coda`

integrating `c++` functionality: `Rcpp`, `RcppArmadillo`

visualization: `ggplot2`, `igraph`

others: `mclust`, `gtools`, `vegan`, `TailRank`




## Usage

To use BiTSC2, please set `R` working directory to SIFA_package after cloning this repository. Make sure you have installed all the dependencies correctly, and then open source code `BiTSC2_app.R` to execute the commands line by line as following.

* For Model Input section, specify random seed `myseed`, and specify the folder `foldername` to save output files (a new folder will be created if it does not exist). For example: 
  ```
  #################################################################
  ###################### MODEL INPUT ##############################
  #################################################################

  myseed <-  1               # set random seed
  foldername <-  "temp_out"          # set output foldername
  dir.create(foldername)  # folder where outputs are saved
  ```

  Then input the total reads matrix and the mutant reads matrix `D` and `X` and squencing depth `psi` and segment information. For the given example data:
  
  ```
  scdata <- readRDS('example_data.RDS')
  D <- scdata$obs.reads$D_drop # total reads, M * N matrix, where row represents locus, column represents cell
  X <- scdata$obs.reads$X_drop # variant reads, M * N matrix. where row represents locus, column represents cell
  #segments <- NULL
  segments <- scdata$segment
  psi <- rep(3,dim(D)[2]) #squencing depth
  ```
If there is genome segment information, it can be used as input information to improve the accuracy of the estimation. If not, assign variable `segment` as `NULL`, that is, use locus specific segments (each gene/ SNV locus as a segment) to update the CNA genotype matrix `L`;
   

* Next, specify Bayesian sampling parameters in `specify_pars.R`. For most of the parameters, default values work just fine. Some of the parameters you can change are:
  ```
  #### maximum number of copy
  Params$max_CN <- 4

  #### maximum number of mutant copies
  Params$max_mut <- 1

  MCMC_par$burnin <- 1000  # burnin sample size
  MCMC_par$Nsamp <- 1000   # number of samples for inference
  MCMC_par$Ntune <- 1000  # number of samples used for adaptive parameter tuning

  Nclone <- c(3:10)  # candidate subclone numbers K
  ```

* Then run the remaining sections one by one:

   > `sampler.R` to perform MCMC sampling, then the samples used for inference are stored in the `temp_out/seed1_K*.Rdata` files;
   
   > `Model_select.R` to make model selection. The corresponding `K` and calculated BIC values are stored in the `temp_out/BIC_model_selection.Rdata`, and the visual graphics of `K` and BIC are displayed in `temp_out/selection.pdf`;
   
   > `Visualization.R` to visualize model sampling results: under different `K`, visualize the estimated subclonal evolutionary tree `T`, CNA genotype matrix `L` and SNV genotype matrix `Z`, and the results are shown in `temp_out/seed1_K*_fit.pdf` files;
   
   > `point_estimate.R` to get the final estimated results from a given posterior sample `.Rdata` file, which are stored in the variable `point_est`.


## Citation
Please cite the BiTSC2 in your publications if it helps your research.
```
@article{chen2020bitsc,
  title={BiTSC\^{} 2: Bayesian inference of Tumor clonal Tree by joint analysis of Single-Cell SNV and CNA data},
  author={Chen, Ziwei and Gong, Fuzhou and Wan, Lin and Ma, Liang},
  journal={bioRxiv},
  year={2020},
  publisher={Cold Spring Harbor Laboratory}
}

```
