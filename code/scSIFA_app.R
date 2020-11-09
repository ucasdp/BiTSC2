
rm(list=ls())
################################################################
############# load packages ##################################
################################################################

setwd('./BiTSC2')

require(tidyr)
require(copynumber)
require(ggplot2)
require(reshape)
require(dplyr)
require(coda)
require(gtools)
require(Rcpp)
require(shape)
require(RcppArmadillo)
require("igraph")
require(mclust)
require(vegan)


########################################################
######## load functions ################################
########################################################

source('./code/assist_fun.R')
source('./code/par_samp.R')
source('./code/tree_samp_fun.R')
source('./code/main_fun.R')
sourceCpp("./code/params_real.cpp")


#################################################################
########## MODEL INPUT ############################################
#################################################################

scdata <- readRDS('example_data.RDS')
# load("example.Rdata")
myseed <-  1               # set random seed
foldername <-  "temp_out"         # set output foldername
dir.create(foldername)  # folder where outputs are saved

D <- scdata$obs.reads$D_drop # total reads, M * N matrix
X <- scdata$obs.reads$X_drop # variant reads, M * N matrix
segments <- NULL
psi <- rep(3,dim(D)[2])


##############################################
######## load parameter file ################
############################################
source('./code/specify_pars.R')
par_disp(Params, MCMC_par)

##############################################
######## sampling ##########################
############################################
source("./code/sampler.R")


##################################################
########## model selection:   ###############
##################################################

source("./code/Model_select.R")
pdf(paste(foldername,"/","selection.pdf",sep=""))
BIC_fun(X,D,foldername)
dev.off()


##################################################
########## visualization   ###############
##################################################
cat("Visualizing sampling results: \n")
source("./code/Visualization.R")
Fit_visual(foldername,X,D)


########################################################
########## get point estimates ###############
########################################################
source("./code/point_estimate.R")
# # specify 
# sample_Rdata <- "seed1_K4.Rdata"
# point_est <- get_point_estimate(foldername,sample_Rdata)


