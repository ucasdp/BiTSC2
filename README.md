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
The steps executed by BiTSC2 are concentrated in the script `BiTSC2_app.R`:

### 1. Input the total reads matrix and the mutant reads matrix `D` and `X` and squencing depth `psi`. 
```
scdata <- readRDS('example_data.RDS')

myseed <-  1               # set random seed
foldername <-  "temp_out"          # set output foldername
dir.create(foldername)  # folder where outputs are saved

D <- scdata$obs.reads$D_drop # total reads, M * N matrix, where row represent locus, column represent cell
X <- scdata$obs.reads$X_drop # variant reads, M * N matrix. where row represent locus, column represent cell
#segments <- NULL
segments <- scdata$segment
psi <- rep(3,dim(D)[2]) #squencing depth
 
 ```
 If there is genome segment information, it can be used as input information to improve the accuracy of the estimation. If not, assign variable `segment` as `NULL`, that is, use locus specific segments (each gene/ SNV locus as a segment) to update the CNA genotype matrix `L`;
2. Initialize the prior parameters;
3. Carry out MCMC sampling, then the samples used for inference are stored in the `temp_out/seed1_K*.Rdata` files;
4. Make model selection. The corresponding `K` and calculated BIC values are stored in the `temp_out/BIC_model_selection.Rdata`, and the visual graphics of `K` and BIC are displayed in `temp_out/selection.pdf`;
5. Visualize model sampling results: under different `K`, visualize the estimated subclonal evolutionary tree `T`, CNA genotype matrix `L` and SNV genotype matrix `Z`, and the results are shown in `temp_out/seed1_K*_fit.pdf` files;
6. According to the model selection, then get the final estimated results, which are stored in the variable `point_est`.


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
