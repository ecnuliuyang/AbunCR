
#################################################################################
###  Abundance estimation when covariates are subject to missing 
###  under discrete time capture-recapture models 
#################################################################################

#################################################################################
###  1. Calculate Ak for a given beta
###  Output: a matrix with m rows and K+1 (or 2K+1) columns
#################################################################################

ak.fun <- function(g.beta, x.obs, K, p){
  Ak <- (1 - g.beta)^K
  if (p == 2) {
    for (k in 1:K) Ak <- cbind( Ak, choose(K, k) * g.beta^k * (1 - g.beta)^(K-k) )
  }else{
    
    for (j in 0:1)
      for (k in 1:K) 
        Ak <- cbind( Ak, choose(K, k) * g.beta^k * (1 - g.beta)^(K-k) * sign(x.obs == j) )
      
  }
  Ak
}


#################################################################################
###  2. minimize  - l_3(nu,beta,lambda) with respect to the Lagrange multiplier 
###  lambda for the given nu and beta
#################################################################################

like3.miss <- function(g.beta, x.obs, n, m, mk, K, nu, p){
  
  el.lam <- function(lam)
    - sum( logg( Ak%*%lam + nu/m) ) - sum( c(nu-n, mk[mk>0]) * log(-lam) )
  
  lambda0 <- c( - (nu - n)/m, - (mk/m)[mk >0] )
  Ak <- ak.fun(g.beta, x.obs, K, p)[, c(TRUE, mk >0)]
  
  nlminb(lambda0 - 1e-5, el.lam, upper = lambda0 - 1e-300)
  
}


#################################################################################
###  3. minimize - l_13(nu,beta) with respect to nu for the given beta 
#################################################################################

like13.miss <- function(g.beta, x.obs, n, m, mk, K, p){
  
  like.nu <- function(nu){
    
    like1 <- sum( log( (nu - n + 1):nu ) ) + (nu - n) * log( (nu - n)/m + 1e-300 )
    like3 <- like3.miss(g.beta, x.obs, n, m, mk, K, nu, p)$objective
    - (like1 + like3)
    
  }
  
  optimize(like.nu, interval=c(n, 100*n), tol=0.01)
  
}


#################################################################################
###  4. Calculate empirical log-likelihood function for the given nu_0  
#################################################################################

el.miss.null <- function(z.mat, x.obs, d.obs, n, m, mk, K, nu0, beta.initial){
  
  like1.null <- sum( log( (nu0 - n + 1):nu0 ) ) + (nu0 - n) * log( (nu0 - n + 1e-300)/m )
  
  like.beta.null <-function(beta){
    
    g.beta <- expit( z.mat%*%beta )
    like2 <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log( 1 - g.beta + 1e-300 ) )
    like3 <- like3.miss(g.beta, x.obs, n, m, mk, K, nu0, length(beta))$objective
    - ( like2 + like3 )
  }
  
  temp <- nlminb(beta.initial, like.beta.null)
  like1.null - temp$objective
  
}



#################################################################################
###  5. Implement empirical likelihood method    
#################################################################################
el.miss <- function(x, z.mat, d, K, beta.init, level){
  
  n <- length(d)
  m <- nrow(z.mat)
  p <- ncol(z.mat)
  d.obs <- d[1:m]; d.mis <- d[(m+1):n]
  x.obs <- x[1:m]; x.mis <- x[(m+1):n]
  
  ##############    definition of mk   ##############
  
  if (p == 2) {
    
    mk <- rep(0, K)
    for (k in 1:K) mk[k] <- sum(d.mis == k)
    
  }else{
    mk <- matrix(0, K, 2)
    for(j in 0:1)
      for(k in 1:K) mk[k, j+1] <- sum( d.mis == k & x.mis == j )
      
      mk <- as.vector(mk)
  }
  
  ##############  EL estimator of beta  ##############
  
  el.beta <- function(beta){
    
    g.beta <- expit( z.mat %*% beta )
    like2.miss <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log( 1 - g.beta + 1e-300 ) )
    temp <- like13.miss(g.beta, x.obs, n, m, mk, K, p)$objective
    - ( like2.miss - temp )
  }
  
  out <- nlminb(beta.init, el.beta)
  like <- - out$objective
  beta.est <- out$par
  g.beta <- expit(z.mat%*%beta.est)
  
  ###  Maximum emiprical log-Likelihood & AIC  ##############
  loglike <- like - sum(log(1:n)) - m * log(m) + sum( mk[mk>0] * log(mk[mk>0]/m) )
  num.par <- 2 + sum(mk>0) + length(beta.est)
  
  ##############  EL estimator of nu  ##############
  out13 <- like13.miss(g.beta, x.obs, n, m, mk, K, p)
  nu.est <- out13$minimum
  
  ############## EL estimator of lambda and alpha  ############
  out1 <- like3.miss(g.beta, x.obs, n, m, mk, K, nu.est, p)
  lam.est <- out1$par
  alpha.est <- - c(nu.est - n, mk[mk>0])/(m*lam.est)
  
  ###########   probability weights of EL  ##############
  prob.est <- 1/(m * (1 + ( ak.fun(g.beta, x.obs, K, p)[, c(TRUE, mk>0)] - 
                              matrix(alpha.est, m, length(alpha.est), byrow = T) )%*%lam.est))
  
  ###########      Test of convergence     ##############
  converge <- ifelse( abs(1 - sum(prob.est) ) < 1e-2 & all(prob.est > 0) & all(prob.est < 1) & 
                        all(ak.fun(g.beta, x.obs, K, p)[, c(TRUE, mk>0)]%*%lam.est + nu.est/m > 0), 0, 1 )
  
  ##############  EL estimator of alpha  ##############
  temp.alpha <- alpha.est[1]
  ind <- 1
  
  if (p == 2) {
    
    for (k in 1:K)
      if ( mk[k] > 0 ) {
        
        ind <- ind + 1
        temp.alpha <- c(temp.alpha, alpha.est[ind])
        
      }else 
        temp.alpha <- c(temp.alpha, choose(K, k) * sum( g.beta^k * (1 - g.beta)^{K - k} * prob.est))
      
  }else{
    
    for (j in 0:1)
      for (k in 1:K)
        if ( mk[K*j + k] > 0 ) {
          
          ind <- ind + 1
          temp.alpha <- c(temp.alpha, alpha.est[ind])
          
        }else
          temp.alpha <- c(temp.alpha, choose(K, k) * sum( g.beta^k * (1 - g.beta)^{K - k} * sign(x.obs == j) * prob.est))
        
  }
  
  alpha.est <- temp.alpha
  
  ###########           ELR based CI      ##############
  rn <- function(nu){
    
    like.full <- like
    like.null <- el.miss.null(z.mat, x.obs, d.obs, n, m, mk, K, nu, beta.est)
    2 * (like.full - like.null) - qchisq(level, 1)
    
  }
  
  ntemp <- nu.est
  while(rn(ntemp) <= 0) ntemp <- ntemp*2
  
  ci.upper <- uniroot(rn, c(nu.est, ntemp))$root
  if (rn(n) <= 0) {
    ci.lower <- n
  }else{
    ci.lower <- uniroot(rn, c(n, nu.est))$root
  }
  
  ci.est <- c(ci.lower, ci.upper)
  
  list(nu.est = nu.est, ci.est = ci.est, alpha.est = alpha.est, beta.est = beta.est, 
       loglike = loglike, num.par = num.par, n = n, m = m, converge = converge)
  
}

#################################################################################
###  6. Calculate empirical log-likelihood under the hypothesis 
###     H0: there is no effect of nonmissing covariates on the capture probability
#################################################################################


loglike.null <- function(cov.bin, cov.miss, cap.times, ind.obs, K){
 
  ### clear the data set
  x <- c(cov.bin[ind.obs == 1], cov.bin[ind.obs == 0])
  d <- c(cap.times[ind.obs == 1], cap.times[ind.obs == 0])
  n <- length(ind.obs)
  m <- sum(ind.obs)
  
  z.mat <- cbind(1, cov.bin[ind.obs == 1], cov.miss[ind.obs == 1])
  rz <- apply(z.mat, 2, max)
  mrz <- matrix(rz, nrow = m, ncol = ncol(z.mat), byrow=T)
  newz <- z.mat/mrz
  
  p <- ncol(newz)
  d.obs <- d[1:m]; d.mis <- d[(m+1):n]
  x.obs <- x[1:m]; x.mis <- x[(m+1):n]
  
  ##############    definition of mk   ##############
  if (p == 2) {
    
    mk <- rep(0, K)
    for (k in 1:K) mk[k] <- sum(d.mis == k)
    
  }else{
    mk <- matrix(0, K, 2)
    for(j in 0:1)
      for(k in 1:K) mk[k, j+1] <- sum( d.mis == k & x.mis == j )
      
      mk <- as.vector(mk)
  }
  
  ##############  EL estimator of beta  ##############
  if (p == 2){
    
    el.beta <-function(beta){
      
      beta <- c(beta, 0)
      g.beta <- expit(newz%*%beta)
      like2.miss <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log( 1 - g.beta + 1e-300 ) )
      tmp <- like13.miss(g.beta, x.obs, n, m, mk, K, p)$objective
      - ( like2.miss - tmp )
    }
    
    out <- optimize(el.beta, c(-1e3, 1e3) )$objective
    
  }else{
    
    el.beta <-function(beta0){
      beta <- rep(0, p)
      beta[1] <- beta0[1]
      beta[p] <- beta0[2]
    
      g.beta <- expit(newz%*%beta)
      like2.miss <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log( 1 - g.beta + 1e-300 ) )
      tmp <- like13.miss(g.beta, x.obs, n, m, mk, K, p)$objective
      - ( like2.miss - tmp )
    }
  
    out <- nlminb(rep(0, 2), el.beta, lower = rep(-1e3, 2), upper = rep(1e3, 2))$objective
 
  }
  
  - out - sum(log(1:n)) - m * log(m) + 
    sum( mk[mk>0] * log(mk[mk>0]/m) )
  
}

