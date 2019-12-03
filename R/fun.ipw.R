#################################################################################
###  The inverse probability weighting and multiple imputation methods
#################################################################################

#################################################################################
###  1. Function to implement the inverse probability weighting method
#################################################################################

ipw.mar <- function(d, K, x = NULL, y, level =0.95, beta.initial = NULL){

  ##############################################################################
  ############   Preparation  ################
  y.mat <- as.matrix(y)
  na.loc <- apply(y.mat, 1, function(y) any(is.na(y)))
  y.obs <- as.matrix(y.mat[!na.loc,])

  d <- as.numeric(d)
  n <- length(d)
  m <- sum(!na.loc)
  d.obs <- d[!na.loc]
  d.mis <- d[na.loc]
  d.all <- c( d.obs, d.mis)

  if ( is.null(x) ) x.mat <- matrix(1, nrow=n)
  else  x.mat <- as.matrix(x)

  x.mis <- as.matrix(x.mat[na.loc,])
  x.obs <- as.matrix(x.mat[!na.loc,])
  x.all <- rbind(x.obs, x.mis)

  ### z.obs = z.mat
  z.obs = cbind(x.obs, y.obs)
  set.seed(321)
  if (is.null(beta.initial)) beta.initial <- runif( ncol(z.obs), -1, 1 )

  z.mat <- as.matrix(z.obs[, beta.initial!=0])
  cov.all <- rbind(z.mat, cbind(x.mis, matrix(0, nrow = n-m, ncol = ncol(y.mat))
                                )[, beta.initial!=0] )
  cov.all <- as.matrix(cov.all)
  ind.obs <- c( rep(1, m), rep(0, n-m) ) ### r


  ##############################################################################
  ############   Main functions  ################
  ### selection probability
  pi.yx.non.para.fun <- function(){

    out <- NULL
    for (i in 1:n) {

      n.obs <- sum( d.obs == d.all[i] &
                      apply( x.obs, 1, function(x) all(x == x.all[i, ]) ) )
      n.all <- sum( d.all == d.all[i] &
                      apply( x.all, 1, function(x) all(x == x.all[i, ]) ) )

      out <- c( out, n.obs/n.all )

    }

    out
  }

  pi.yx <- pi.yx.non.para.fun()

  pi.yx.non.na <- 0
  for (i in 1:n) pi.yx.non.na[i] <- ifelse ( pi.yx[i] == 0, 1e-300, pi.yx[i] )

  z.max <- apply(z.mat, 2, max)
  z.matm <- (z.mat/matrix(z.max, nrow = nrow(z.mat), ncol = ncol(z.mat), byrow = T))

  like.beta.ipw <- function (beta) {
    beta <- as.matrix(beta)
    g.beta <- as.numeric( plogis( z.matm %*% beta ) )
    phi.beta <- 1 - (1 - g.beta)^K
    - sum( (d.obs * log(g.beta + 1e-300) + (K - d.obs)*log(1 - g.beta + 1e-300) -
              log(phi.beta + 1e-300))/pi.yx[1:m] )
  }

  out <- nlminb(beta.initial[beta.initial!=0],
                like.beta.ipw)
  beta.est <- as.matrix(out$par/z.max)
  gi <- as.numeric( plogis(cov.all%*%beta.est) )
  p.star <- 1 - ( 1 - gi )^K
  nu.est <- sum( ind.obs/( pi.yx.non.na * p.star ) )


  ### Variance estimator of beta
  psi.i <- NULL
  p <- nrow(beta.est)
  psi.i.star <- matrix(0, n, p)
  kap.star <- rep(0, n)
  for (i in 1:n)
    psi.i <- rbind( psi.i, (d.all[i] - K*gi[i]/p.star[i]) * matrix(cov.all[i,],nrow = 1) )
  for (i in 1:n) {

    m.00 <- 0; m.0 <- 0; m.1 <- 0; m.2 <- rep(0, p)

    for ( j in 1:n ) {

      if ( d.all[j] == d.all[i] ){

        m.00 <- m.00 + 1
        m.0 <- m.0 + ind.obs[j] / pi.yx.non.na[j]
        m.1 <- m.1 + ind.obs[j] / ( pi.yx.non.na[j] * p.star[j] )
        m.2 <- m.2 + ind.obs[j] * psi.i[j,] / pi.yx.non.na[j]

      }

    }

    if(m.0 > 0){
      kap.star[i] <- m.1/m.0
      psi.i.star[i,] <- m.2/m.00
    }
  }

  Uwi <- psi.i * ind.obs / pi.yx.non.na + psi.i.star * ( 1 - ind.obs / pi.yx.non.na )

  Gw <- 0
  for (i in 1:m)  Gw <-  Gw - ( K^2 * gi^2 * (1 - p.star) - K * gi * (1-gi) * p.star )[i] *
    as.matrix(cov.all[i,]) %*% matrix(cov.all[i,], nrow=1) / ( p.star[i]^2 * pi.yx[i] )

  var.beta <- solve(Gw) %*% (
    t(Uwi) %*% Uwi ) %*% t(solve(Gw))

  ### Variance estimator of nu
  A.hat <- t( - t(cov.all) %*% ( ind.obs * K * gi * (1 - p.star) / (pi.yx.non.na * p.star^2) ) ) %*% t( solve(Gw) )

  out <- ind.obs / (pi.yx.non.na * p.star) - kap.star * (ind.obs - pi.yx.non.na) / pi.yx.non.na +
    ( ind.obs * psi.i / pi.yx.non.na - (ind.obs - pi.yx.non.na) * psi.i.star) %*% t(A.hat)
  sigma2.nu <- as.numeric( sum(out^2) - nu.est )


  rt <- list(n.big = nu.est, n.big.se = sqrt(sigma2.nu),
             n.big.ci = c( nu.est - sqrt(sigma2.nu)*qnorm(0.5+level/2),
                           nu.est + sqrt(sigma2.nu)*qnorm(0.5+level/2) ),
             beta = as.numeric(beta.est), beta.se =  sqrt(diag(var.beta) ),
             level = level)
  class(rt) <- 'IPW'
  return(rt)
}


print.IPW <- function (obj) {

  cat ("\nInverse probability weighting method: \n")
  cat ("\nCoefficients:\n")
  beta <- matrix(round(obj$beta, 2),nrow=1)
  coefficients <- rbind( beta, round(obj$beta.se[1:ncol(beta)], 2) )
  rownames(coefficients) <- c('Estimate','Std. Error')
  colnames(coefficients) <- paste0('beta', 1:ncol(coefficients))

  print(coefficients)

  cat ("\nAbundance:\n")
  cat ("Estimate:", round(obj$n.big))
  cat ("\nStd. Error:", round(obj$n.big.se))

  cat ("\n")
  cat (paste0(100*obj$level,"%"),
       "CI: [",
       paste(round(obj$n.big.ci), collapse=', '), "]\n")


}

