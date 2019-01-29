#################################################################################
###               Predefined functions      
#################################################################################

#################################################################################
###  1. Function to approximate log(t)          
#################################################################################

logg <- function(t){
  
  eps <- 1e-5
  ts <- t*(t>eps)+1*(t<=eps)
  log(ts) * (t>eps) + (log(eps) - 1.5 + 2*t/eps - (t/eps)^2*0.5) * (t<=eps)
  
}


#################################################################################
###  2. For the given x, calculate exp(x)/(1+exp(x) )                     
#################################################################################

expit <- function(x)  plogis(x, 0, 1)


#################################################################################
###   Main function for capture recapture model when covariates are subject to missing
###         Implement the empirical likelihoods & IPW & MI methods
#################################################################################


CRMiss <- function(cov.bin = NULL, cov.miss, cap.times, ind.obs, K, level = 0.95, M = 100){
  
  ### clear the data set
  x <- c(cov.bin[ind.obs == 1], cov.bin[ind.obs == 0])
  d <- c(cap.times[ind.obs == 1], cap.times[ind.obs == 0])
  n <- length(ind.obs)
  m <- sum(ind.obs)
  
  z.mat <- cbind(1, cov.bin[ind.obs == 1], cov.miss[ind.obs == 1])
  rz <- apply(z.mat, 2, max)
  mrz <- matrix(rz, nrow = m, ncol = ncol(z.mat), byrow=T)
  newz <- z.mat/mrz
  
  ###   The inverse probability weighting method
  
  out.ipw <- ipw.miss(cov.bin[ind.obs == 0], newz, cap.times[ind.obs == 1], cap.times[ind.obs == 0], K)
  
  ###   The multiple imputation method
  
  out.mi <- im2.miss(cov.bin[ind.obs == 0], newz, cap.times[ind.obs == 1], cap.times[ind.obs == 0], K, M)
  
  ###   The naive empirical likelihood method
  
  out.el.cc <- el.cc(newz, cap.times[ind.obs == 1], K, beta.init = out.ipw$beta.est, level)
  
  ###   The proposed empirical likelihood method
  
  out.el <- el.miss(x, newz, d, K, beta.init = out.ipw$beta.est, level)
  
  ### Restore the estimates of beta
  
  out.ipw$beta.est <- out.ipw$beta.est/rz
  out.mi$beta.est <- out.mi$beta.est/rz
  out.el.cc$beta.est <- out.el.cc$beta.est/rz
  out.el$beta.est <- out.el$beta.est/rz
  
  ### confidence intervals of nu
  temp.ipw <- qnorm(0.5 + level/2) * sqrt(out.ipw$sigma2.nu)
  temp.mi <- qnorm(0.5 + level/2) * sqrt(out.mi$sigma2.nu)
  
  out.ipw$ci.est <- c( out.ipw$nu.est - temp.ipw, out.ipw$nu.est + temp.ipw ) 
  out.mi$ci.est <- c( out.mi$nu.est - temp.mi, out.mi$nu.est + temp.mi ) 
  
  ### estimates of alpha
  if ( is.null(cov.bin) ) {
  
       temp.alpha <- out.el$alpha.est[2:(K+1)]
       names(temp.alpha) <- paste0('gamma', 1:K)
       
  }else{
  
       temp.alpha <- matrix(out.el$alpha.est[2:(2*K+1)], 2, K, byrow = T)
       colnames(temp.alpha) <- paste0('k=', 1:K)
       rownames(temp.alpha) <- c('j=0', 'j=1')
       
  }
  
  ### Output
  list(num.occassion = K, num.cap = n, num.obs = m, 
       nu.el = out.el$nu.est, ci.el = out.el$ci.est, 
       beta.el = out.el$beta.est, gamma0.el = out.el$alpha.est[1], gammajk.el = temp.alpha,
       loglike.el = out.el$loglike, num.para = out.el$num.par, converge.el = out.el$converge,
       nu.el.cc = out.el.cc$nu.est, ci.el.cc = out.el.cc$ci.est, beta.el.cc = out.el.cc$beta.est, 
       nu.ipw = out.ipw$nu.est, sigma2.nu.ipw = out.ipw$sigma2.nu, ci.ipw = out.ipw$ci.est, beta.ipw = out.ipw$beta.est, 
       nu.mi = out.mi$nu.est, sigma2.nu.mi = out.mi$sigma2.nu, ci.mi = out.mi$ci.est, beta.mi = out.mi$beta.est)
  
}

#################################################################################
###  summary statatistics of (beta, nu)
#################################################################################
nu.summary <- function(outcomes){
  
  out <- matrix(0, 4, 3)
  rownames(out) <- c('EL', 'IPW', 'MI', 'EL-CC')
  colnames(out) <- c('est', 'lower', 'upper')
  
  out[1,1] <- outcomes$nu.el
  out[2,1] <- outcomes$nu.ipw
  out[3,1] <- outcomes$nu.mi
  out[4,1] <- outcomes$nu.el.cc
  out[1,2:3] <- outcomes$ci.el
  out[2,2:3] <- outcomes$ci.ipw
  out[3,2:3] <- outcomes$ci.mi
  out[4,2:3] <- outcomes$ci.el.cc
  
  round(out)
  
}

beta.summary <- function(outcomes){
  
  p <- length(outcomes$beta.ipw)
  out <- matrix(0, 4, p)
  rownames(out) <- c('EL', 'IPW', 'MI', 'EL-CC')
  colnames(out) <- paste0('beta', 1:p)
  
  out[1,] <- outcomes$beta.el
  out[2,] <- outcomes$beta.ipw
  out[3,] <- outcomes$beta.mi
  out[4,] <- outcomes$beta.el.cc
  
  round(out, 4)
  
}




