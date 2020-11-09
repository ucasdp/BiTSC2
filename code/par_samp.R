
#####################################################
########## sampling of pi #########################
#######################################################

samp_pi <- function(samp, Params, temper)
{
  L <- samp$L
  temp <- apply(L,1,function(x){all(x==2)})
  n <- sum(temp)
  m <- Params$M-n
  rbeta(1,(n+Params$alpha+temper-1)/temper, (m+Params$beta+temper-1)/temper)
}

#####################################################
########## sampling of rho #########################
#######################################################

samp_rho <- function(samp, Params, D, temper)
{
  rho0 <- samp$rho
  rho <- rho0
  LL <- samp$L[,samp$C]

  psi <- Params$psi
  rho1 <- runif(1,min=0,max=1)
  p0 <- sum(log_prob_D(D, LL, psi, rho0))/temper
  p1 <- sum(log_prob_D(D, LL, psi, rho1))/temper
  alpha <- exp(p1-p0)
  u <- runif(1,min=0,max=1)
  if(u<min(1,alpha)){rho <- rho1}
  rho
}


##################################################################
########## tune theta parameter ###############################
##################################################################
adaptive_tune <- function(S,low=0.4,high=0.6)
{
  mcmcS <-  mcmc(S)
  rej_rate <- rejectionRate(mcmcS)
  rej_rate <- mean(rej_rate)
  
  if(rej_rate>high)
  {out <- 1} else if (rej_rate<low)
  {out <- -1} else {out <- 0}
  out2 <- list(move=out,rate=rej_rate)
  out2
}


