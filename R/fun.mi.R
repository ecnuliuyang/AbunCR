#################################################################################
###  The inverse probability weighting and multiple imputation methods
#################################################################################

#################################################################################
###  2. Function to implement the multiple imputation method
#################################################################################

mi2.mar <- function(d, K, x = NULL, y, level = 0.95, M = 100, seed = 321, beta.initial = NULL){

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

  x.mis <- matrix(x.mat[na.loc,], ncol = ncol(x.mat))
  x.obs <- matrix(x.mat[!na.loc,], ncol = ncol(x.mat))
  x.all <- rbind(x.obs, x.mis)


  ### z.obs = z.mat
  z.obs = cbind(x.obs, y.obs)
  set.seed(321)
  if (is.null(beta.initial)) beta.initial <- runif( ncol(z.obs), -1, 1 )

  p <- ncol(z.obs)
  z.mat <- as.matrix(z.obs[, beta.initial!=0])
  cov.all <- rbind(z.mat, cbind(x.mis, rep(0, n-m))
                   [, beta.initial!=0] )
  cov.all <- as.matrix(cov.all)

  ind.obs <- c( rep(1, m), rep(0, n-m) ) ### r

  data <- cbind(d.all, x.all, c(y.obs, rep(0, n-m)), c(rep(1, m), rep(0, n-m)) )
  data <- as.matrix(data)
  y.max <- max(y.obs)


  ##############################################################################
  ############   Main functions  ################

  MI.data.fun <- function(){

    mi.data <- matrix(0, n, M+1)
    colnames(mi.data) <- c( paste0('Z',1:M), 'R' )
    mi.data[, M+1] <- data[, ncol(data)]
    for (i in 1:n)
      if (data[i, ncol(data)] == 1) {
        mi.data[i, 1:M] <- rep(y.obs[i], M)
      }else{

        xx <- y.obs[ d.obs == d.all[i] & apply( x.obs, 1, function(x) all(x == x.all[i,]) ) ]
        if (length(xx) != 0) mi.data[i, 1:M] <- sample(xx, M, replace = T)
        else mi.data[i,1:M] <- rep( mean(y.obs), M)

      }


    mi.data
  }

  like.con.mi <- function(beta){

    beta <- as.matrix(beta)
    like.i <- function(d, p.i, p.i.star)
      d*log(p.i + 1e-300) + (K - d)*log(1 - p.i + 1e-300) - log(p.i.star + 1e-300)

    out <- 0
    for (mm in 1:M) {
      dat.m[[mm]][,ncol(dat.m[[mm]])] <- dat.m[[mm]][,ncol(dat.m[[mm]])]/y.max
      p.i <- plogis( dat.m[[mm]]%*%beta )
      p.i.star <- 1 - (1 - p.i)^K
      out <- out + sum( like.i(d.all, p.i, p.i.star) )

    }
    - out/M
  }

  df.2 <- function(beta0){

    out <- 0
    for(mm in 1:M){
      xx <- dat.m[[mm]]
      p.i <- plogis(xx%*%beta0)
      p.i.star <- 1 - (1-p.i)^K
      item <- - as.numeric( (K * p.i * (1-p.i) * p.i.star - K^2 * p.i^2 * (1-p.i.star) )/(p.i.star^2) )
      out <- out + t(xx) %*% diag(item) %*% xx
    }
    out/M
  }

  N.est.2 <- function(beta){

    p.i.star.mi <- 0

    for (mm in 1:M) {
      p.i <- plogis( dat.m[[mm]]%*%beta )
      p.i.star <- 1 - (1 - p.i)^K
      p.i.star.mi <- p.i.star.mi + 1/p.i.star
    }
    sum(1/(M/p.i.star.mi))
  }

  var.est.2 <- function(beta){

    # beta <- beta.est
    GG <- - df.2(beta)

    u.nu.i <- 0; u.ni <- 0; v.1 <- 0; v.3 <- 0; nu.est.mi <- NULL
    for (mm in 1:M) {

      xx <- dat.m[[mm]]
      p.i <- plogis(xx%*%beta)
      p.i.star <- 1 - (1 - p.i)^K
      u.nu.i <- u.nu.i + t(xx) %*% diag(as.numeric( (d.all - K*p.i/p.i.star)^2 )) %*% xx
      u.ni <- u.ni + t(xx) %*% (d.all - K*p.i/p.i.star) %*% t( t(xx) %*% (d.all - K*p.i/p.i.star) )
      v.1 <- v.1 + sum( (1-p.i.star)/p.i.star^2 )
      v.3 <- v.3 - t(xx) %*% ( K*(1 - p.i.star)*p.i/(p.i.star^2) )
      nu.est.mi <- c( nu.est.mi, sum(1/p.i.star) )

    }

    vv <- solve(GG) %*% (u.nu.i/M + (1+1/M) * u.ni / (M-1)) %*% solve(GG)

    c( diag(vv), v.1/M + (1+1/M) * sum( (nu.est.mi - est.mi[length(est.mi)])^2 )/(M-1) + t(v.3/M) %*% vv %*% (v.3/M) )

  }

  set.seed(seed)
  dat.m <- list()

  mi.data <- MI.data.fun()
  for (mm in 1:M)
    dat.m[[mm]] <- cbind(x.all, mi.data[,mm])[, beta.initial!=0]

  out <- nlminb(beta.initial[beta.initial!=0], like.con.mi)
  beta.est <- out$par
  beta.est[length(beta.est)] <- beta.est[length(beta.est)]/y.max

  est.mi <- c(beta.est, N.est.2(beta.est))
  est.var <- var.est.2( beta.est )

  rt <- list(n.big = est.mi[length(est.mi)],
             n.big.se = sqrt(est.var[length(est.var)]),
             n.big.ci = c( est.mi[length(est.mi)] - sqrt(est.var[length(est.var)])*qnorm(0.5+level/2),
                           est.mi[length(est.mi)] + sqrt(est.var[length(est.var)])*qnorm(0.5+level/2) ),
             beta = est.mi[1:(length(est.mi) - 1)],
             beta.se =  sqrt(est.var[1:(length(est.var)-1)]),
             level = level)
  class(rt) <- 'MI'
  return(rt)

}


print.MI <- function (obj) {

  cat ("\nMultiple Imputation method: \n")
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





