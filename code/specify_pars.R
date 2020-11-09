################################################################################
#
# this file is used to specify model parameters
#
########################################################################################


########################################################
############## MODEL PARAMETERS ############################
########################################################

set.seed(myseed)

MCMC_par <- Params <- list()

#### number of cells
Params$N <-  ncol(D)

#### number of genes
Params$M <- nrow(D)

#### maximum number of copy
Params$max_CN <- 4

#### maximum number of mutant copies
Params$max_mut <- 2

##### hyper parameter of phi: theta
Params$r=1.5

##### hyper parameter of pi
Params$alpha <- 10^4
Params$beta <- 1

##### hyper parameter in prior of Z: prob of gaining one mutant allele
Params$zeta <- 10^-2

##### read depth parameter
Params$psi <- psi



########################################################
############## get segments ############################
########################################################

if(is.null(segments)){
  cat("Genome segment information not found. Generating segments with bin size 1\n")
  segments <- cbind(c(1:dim(D)[1]),c(1:dim(D)[1]))
}

Params$segments <- segments


########################################################
############  specify MCMC parameters ##################
########################################################

MCMC_par$burnin <- 500  # burnin sample size
MCMC_par$Nsamp <- 500   # number of samples for inference
MCMC_par$Ntune <- 500  # number of samples used for adaptive parameter tuning



MCMC_par$swap_interval <- 30 # make Matroplis Hastings move in every how many samples
MCMC_par$Nchain <- 5 # number of paralel chains
MCMC_par$delta_T <- 0.35

Temperature <- seq(1,by=MCMC_par$delta_T,length.out = MCMC_par$Nchain )  # temperatures

# ######################################################
# ############## adaptive tuning parameter ##############
# ######################################################
 
adapt_par <- list(theta_tune=6)
adapt_par <- lapply(c(1:MCMC_par$Nchain),function(x){adapt_par})


######################################################
##########  sequence of candidate Ks ##########
######################################################

# candidate subclone numbers K
# need K >= 2
Nclone <- c(3:7) 

