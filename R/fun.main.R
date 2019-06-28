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

  ### tidy the estimates of beta

  p <- length(out.ipw$beta.est)
  out <- matrix(0, 4, p)
  rownames(out) <- c('MELE', 'IPW', 'MI', 'CC')
  colnames(out) <- paste0('beta', 1:p)

  out[1,] <- out.el$beta.est
  out[2,] <- out.ipw$beta.est
  out[3,] <- out.mi$beta.est
  out[4,] <- out.el.cc$beta.est

  out.beta <- round(out, 4)

  ### confidence intervals of nu

  temp.ipw <- qnorm(0.5 + level/2) * sqrt(out.ipw$sigma2.nu)
  temp.mi <- qnorm(0.5 + level/2) * sqrt(out.mi$sigma2.nu)

  out.ipw$ci.est <- c( out.ipw$nu.est - temp.ipw, out.ipw$nu.est + temp.ipw )
  out.mi$ci.est <- c( out.mi$nu.est - temp.mi, out.mi$nu.est + temp.mi )

  ### tidy the estimates of nu
  out <- matrix(0, 4, 3)
  rownames(out) <- c('MELE', 'IPW', 'MI', 'CC')
  colnames(out) <- c('est', 'lower', 'upper')

  out[1,1] <- out.el$nu.est
  out[2,1] <- out.ipw$nu.est
  out[3,1] <- out.mi$nu.est
  out[4,1] <- out.el.cc$nu.est
  out[1,2:3] <- out.el$ci.est
  out[2,2:3] <- out.ipw$ci.est
  out[3,2:3] <- out.mi$ci.est
  out[4,2:3] <- out.el.cc$ci.est

  out.nu <- round(out)

  ### tidy the maximum empirical likelihood estimate of alpha
  if ( is.null(cov.bin) ) {

       temp.alpha <- out.el$alpha.est[2:(K+1)]
       names(temp.alpha) <- paste0('gamma', 1:K)

  }else{

       temp.alpha <- matrix(out.el$alpha.est[2:(2*K+1)], 2, K, byrow = T)
       colnames(temp.alpha) <- paste0('k=', 1:K)
       rownames(temp.alpha) <- c('j=0', 'j=1')

  }

  alpha.el <- list(gamma0.el = out.el$alpha.est[1],
                   gammajk.el = temp.alpha)

  ### Output
  list(nu.est = out.nu, beta.est = out.beta, alpha.el = alpha.el,
       loglike.el = out.el$loglike, converge = out.el$converge)

}

