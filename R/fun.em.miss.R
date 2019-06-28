
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
  
  ###########      Variance estimates      ##############
  if(p == 2){
    
    if(sum(mk!=0) == K) se <- sqrt(sig2.mele.g1(z.mat, d.obs, d.mis, mk, nu.est, beta.est, lam.est, alpha.est, prob.est))[1:(p+1)]
    if(sum(mk!=0) < K) se <- sqrt(sig2.mele.s1(z.mat, d.obs, d.mis, mk, nu.est, beta.est, lam.est, alpha.est, prob.est))[1:(p+1)]
    
  }else{
    
    if(sum(mk!=0) == K*2) se <- sqrt(sig2.mele.g2(z.mat, d.obs, d.mis, x.obs, x.mis, mk, nu.est, beta.est, lam.est, alpha.est, prob.est))[1:(p+1)]
    if(sum(mk!=0) < K*2) se <- sqrt(sig2.mele.s2(z.mat, d.obs, d.mis, x.obs, x.mis, mk, nu.est, beta.est, lam.est, alpha.est, prob.est))[1:(p+1)]
    
  }
  
  
  
  
  
  list(nu.est = nu.est, ci.est = ci.est, alpha.est = alpha.est, beta.est = beta.est, se.nu = se[1], se.beta = se[-1],
       loglike = loglike, num.par = num.par, n = n, m = m, converge = converge)
  
}


#################################################################################
###  6. Functions to calculate the variance estimates of MELE
#################################################################################

sig2.mele.s1 <- function(z.mat, d.obs, d.mis, mk, nu, beta, lam, alp, pi){

  g.beta <- plogis(z.mat%*%beta)
  Ak <- ak.fun(g.beta, x.obs = NULL, K, p)
  ind_pos <- mk > 0
  Ak_pos <- Ak[, c(TRUE, ind_pos)]
  
  hk <- rep(1, K)
  for (k in 1:K){
    if (ind_pos[k]) hk[k] <- sum(d.obs == k)/sum( c(d.obs, d.mis) == k )
  } 
  
  lam00 <- sum(hk*alp[-1])
  hk.pos <- hk[ind_pos]
  H1s <- as.matrix( c(-1, hk.pos - 1) )
  H2s <- diag( c(1/alp[1], (1-hk.pos)/alp[-1][ind_pos]) )
  
  pis <- 0
  for (k in 1:K) pis <- pis + choose(K, k) * 
    g.beta^k * (1-g.beta)^{K-k} * hk[k]
  
  V11 <- 1 - 1/alp[1]
  V31 <- as.matrix( 1/alp[1]*c(1, rep(0, sum(ind_pos))) )
  V33 <- - H2s + 1/nu*sum(1/pis^2)*H1s%*%t(H1s)
  
  V22 <- 0
  V23 <- 0
  V43 <- 0
  V24 <- 0
  V44 <- 0
  
  
  Uk <- (1 - g.beta)^K - alp[1]
  for(k in which(ind_pos)){
    Uk <- cbind( Uk, choose(K, k) * g.beta^k * (1 - g.beta)^(K-k) - alp[k+1] )
  }
  
  for (i in 1:m) {
    pi.dot <- 0
    pi.ddot <- 0
    for (k in 1:K){
      pi.dot <- pi.dot + choose(K, k) * 
        g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[k] *
        (k - K*g.beta[i]) * as.matrix(z.mat[i,]) 
      
      pi.ddot <- pi.ddot + choose(K, k) * 
        g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[k]*
        ( (k - K*g.beta[i])^2 - K*g.beta[i]*(1-g.beta[i]) ) * 
        as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) )
    }
    
    
    V22 <- V22 + ( pi.dot%*%t(pi.dot)/pis[i] - pi.ddot - 
                     K* g.beta[i] * (1-g.beta[i]) * pis[i] *
                     as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) ) )/pis[i]
    V23 <- V23 - pi.dot/pis[i]^2
    V43 <- V43 + Uk[i,]/pis[i]^2
    
    # i=1
    Uk.dot <- (1 - g.beta[i])^K * (0 - K*g.beta[i]) * z.mat[i, ]
    for(k in which(ind_pos)) 
      Uk.dot <- cbind( Uk.dot, choose(K, k) * g.beta[i]^k * (1 - g.beta[i])^(K-k) * 
                         (k - K*g.beta[i]) * z.mat[i, ] )
    
    V24 <- V24 - Uk.dot/pis[i] + pi.dot%*%t(as.matrix(Uk[i,]))/pis[i]^2
    V44 <- V44 + as.matrix(Uk[i,])%*%t(as.matrix(Uk[i,]))/pis[i]^2
  }
  V22 <- V22/nu
  V23 <- as.matrix(V23/nu)%*%t(H1s)
  V43 <- lam00*( diag(1,sum(ind_pos)+1) - as.matrix(V43/nu)%*%t(H1s) )
  V24 <- lam00*V24/nu
  V44 <- lam00^2*V44/nu
  
  
  inv.V44 <- solve(V44)
  # inv.V44 <- ginv(V44)
  
  W22 <- - V22 + V24%*%inv.V44%*%t(V24)
  W23 <- - V23 + V24%*%inv.V44%*%V43
  W33 <- - V33 + t(V43)%*%inv.V44%*%V43
  Ws <- rbind( cbind( -V11, matrix(0, 1, p), - t(V31) ),
               cbind( matrix(0, p, 1), W22, W23 ),
               cbind( -V31, t(W23), W33))
  
  ### variance estimates of MELE
  inv.Ws <- solve(Ws)
  # inv.Ws <- ginv(Ws)
  diag(inv.Ws)*c(nu, rep(1/nu,nrow(Ws)-1))
}

sig2.mele.g1 <- function(z.mat, d.obs, d.mis, mk, nu, beta, lam, alp, prob.est){
  
  alp = alpha.est
  hk <- NULL
  for (k in 1:K) hk <- c( hk, sum(d.obs == k)/sum( c(d.obs, d.mis) == k) )
  
  lam00 <- sum(hk*alp[-1])
  H1 <- as.matrix(hk)
  H2 <- diag( c((1-hk)/alp[-1]) )
  
  g.beta <- plogis(z.mat%*%beta)
  m <- nrow(z.mat)
  pi <- 0
  for (k in 1:K) pi <- pi + choose(K, k) * 
    g.beta^k * (1-g.beta)^{K-k} * hk[k]
  
  V11 <- 1 - 1/alp[1]
  V31 <- as.matrix( -1/alp[1]*rep(1, K) )
  V33 <- -1/alp[1]*rep(1, K)%*%t(rep(1, K)) - H2 + 1/nu*sum(1/pi^2)*H1%*%t(H1)
  
  V22 <- 0
  V23 <- 0
  V43 <- 0
  V24 <- 0
  V44 <- 0
  
  Uk <- NULL
  for(k in 1:K){
    Uk <- cbind( Uk, choose(K, k) * g.beta^k * (1 - g.beta)^(K-k) - alp[k+1] )
  }
  
  for (i in 1:m) {
    pi.dot <- 0
    pi.ddot <- 0
    for (k in 1:K){
      pi.dot <- pi.dot + choose(K, k) * 
        g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[k] *
        (k - K*g.beta[i]) * as.matrix(z.mat[i,]) 
      
      pi.ddot <- pi.ddot + choose(K, k) * 
        g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[k]*
        ( (k - K*g.beta[i])^2 - K*g.beta[i]*(1-g.beta[i]) ) * 
        as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) )
    }
    
    
    V22 <- V22 + ( pi.dot%*%t(pi.dot)/pi[i] - pi.ddot - 
                     K* g.beta[i] * (1-g.beta[i]) * pi[i] *
                     as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) ) )/pi[i]
    V23 <- V23 - pi.dot/pi[i]^2
    V43 <- V43 + Uk[i,]/pi[i]^2
    
    Uk.dot <- NULL
    for(k in 1:K) Uk.dot <- cbind( Uk.dot, choose(K, k) * g.beta[i]^k * (1 - g.beta[i])^(K-k) * 
                                     (k - K*g.beta[i]) * z.mat[i, ] )
    
    V24 <- V24 - Uk.dot/pi[i] + pi.dot%*%t(as.matrix(Uk[i,]))/pi[i]^2
    V44 <- V44 + as.matrix(Uk[i,])%*%t(as.matrix(Uk[i,]))/pi[i]^2
  }
  V22 <- V22/nu
  V23 <- as.matrix(V23/nu)%*%t(H1)
  V43 <- lam00*( diag(1,K) - as.matrix(V43/nu)%*%t(H1) )
  V24 <- lam00*V24/nu
  V44 <- lam00^2*V44/nu
  
  inv.V44 <- solve(V44)
  # inv.V44 <- ginv(V44)
  
  W22 <- - V22 + V24%*%inv.V44%*%t(V24)
  W23 <- - V23 + V24%*%inv.V44%*%V43
  W33 <- - V33 + t(V43)%*%inv.V44%*%V43
  W <- rbind( cbind( -V11, matrix(0, 1, p), - t(V31) ),
              cbind( matrix(0, p, 1), W22, W23 ),
              cbind( -V31, t(W23), W33))
  
  ### variance estimates of MELE
  inv.W <- solve(W)
  # inv.W <- ginv(W)
  
  diag(inv.W)*c(nu, rep(1/nu,nrow(W) - 1))
}

sig2.mele.s2 <- function(z.mat, d.obs, d.mis, x.obs, x.mis, mk, nu.est, beta.est, lam.est, alpha.est, prob.est){
  
  # pars <- est_ours
  nu <- nu.est
  beta <- beta.est
  p <- length(beta)
  g.beta <- plogis(z.mat%*%beta)
  m <- nrow(z.mat)
  
  lam <- lam.est
  alp <- alpha.est
  Ak <- ak.fun(g.beta, x.obs, K, p) 
  ind_pos <- mk > 0
  Ak_pos <- Ak[, c(TRUE, ind_pos)]
  pi <- prob.est
  
  
  hk <- NULL
  for (j in 0:1)
    for (k in 1:K) 
      hk <- c( hk, sum(d.obs == k&x.obs == j)/sum(c(d.obs, d.mis) == k&c(x.obs, x.mis) == j) )
  
  hk[is.na(hk)] <- 1
  
  lam00 <- sum(hk*alp[-1])
  hk.pos <- hk[ind_pos]
  H1s <- as.matrix( c(-1, hk.pos - 1) )
  H2s <- diag( c(1/alp[1], (1-hk.pos)/alp[-1][ind_pos]) )
  
  pis <- 0
  for (j in 0:1)
    for (k in 1:K) 
      pis <- pis + choose(K, k) * g.beta^k * (1-g.beta)^{K-k} * 
    (x.obs==j) * hk[j*K+k]
  
  V11 <- 1 - 1/alp[1]
  V31 <- as.matrix( 1/alp[1]*c(1, rep(0, sum(ind_pos))) )
  V33 <- - H2s + 1/nu*sum(1/pis^2)*H1s%*%t(H1s)
  
  V22 <- 0
  V23 <- 0
  V43 <- 0
  V24 <- 0
  V44 <- 0
  
  Uk <- (Ak - matrix(alp,nrow(Ak),ncol(Ak), 
                     byrow = TRUE))[,c(TRUE, ind_pos)]
  
  for (i in 1:m) {
    pi.dot <- 0
    pi.ddot <- 0
    for(j in 0:1)
      for (k in 1:K){
        pi.dot <- pi.dot + choose(K, k) * 
          g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[j*K+k] *(x.obs[i]==j)*
          (k - K*g.beta[i]) * as.matrix(z.mat[i,]) 
        
        pi.ddot <- pi.ddot + choose(K, k) * 
          g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[j*K+k] *(x.obs[i]==j)*
          ( (k - K*g.beta[i])^2 - K*g.beta[i]*(1-g.beta[i]) ) * 
          as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) )
      }
    
    V22 <- V22 + ( pi.dot%*%t(pi.dot)/pis[i] - pi.ddot - 
                     K* g.beta[i] * (1-g.beta[i]) * pis[i] *
                     as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) ) )/pis[i]
    V23 <- V23 - pi.dot/pis[i]^2
    V43 <- V43 + Uk[i,]/pis[i]^2
    
    # i=1
    Uk.dot <- (1 - g.beta[i])^K * (0 - K*g.beta[i]) * z.mat[i, ]
    for (j in 0:1)
      for(k in 1:K) Uk.dot <- cbind( Uk.dot, choose(K, k) * g.beta[i]^k * (1 - g.beta[i])^(K-k) * (x.obs[i]==j)*
                                       (k - K*g.beta[i]) * z.mat[i, ] )
    
    Uk.dot <- Uk.dot[, c(TRUE, ind_pos)]
    
    V24 <- V24 - Uk.dot/pis[i] + pi.dot%*%t(as.matrix(Uk[i,]))/pis[i]^2
    V44 <- V44 + as.matrix(Uk[i,])%*%t(as.matrix(Uk[i,]))/pis[i]^2
  }
  V22 <- V22/nu
  V23 <- as.matrix(V23/nu)%*%t(H1s)
  V43 <- lam00*( diag(1,sum(ind_pos)+1) - as.matrix(V43/nu)%*%t(H1s) )
  V24 <- lam00*V24/nu
  V44 <- lam00^2*V44/nu
  
  
  inv.V44 <- solve(V44)
  # inv.V44 <- ginv(V44)
  
  W22 <- - V22 + V24%*%inv.V44%*%t(V24)
  W23 <- - V23 + V24%*%inv.V44%*%V43
  W33 <- - V33 + t(V43)%*%inv.V44%*%V43
  Ws <- rbind( cbind( -V11, matrix(0, 1, p), - t(V31) ),
               cbind( matrix(0, p, 1), W22, W23 ),
               cbind( -V31, t(W23), W33))
  
  ### variance estimates of MELE
  inv.Ws <- solve(Ws)
  # inv.Ws <- ginv(Ws)
  diag(inv.Ws)*c(nu, rep(1/nu,nrow(Ws)-1))
}

sig2.mele.g2 <- function(z.mat, d.obs, d.mis, x.obs, x.mis, mk, nu.est, beta.est, lam.est, alpha.est, prob.est){
  # pars <- est_ours
  nu <- nu.est
  beta <- beta.est
  g.beta <- plogis(z.mat%*%beta)
  p <- length(beta)
  m <- nrow(z.mat)
  
  lam <- lam.est
  alp <- alpha.est
  
  hk <- NULL
  for (j in 0:1)
    for (k in 1:K) 
      hk <- c( hk, sum(d.obs == k&x.obs == j)/sum(c(d.obs, d.mis) == k&c(x.obs, x.mis) == j) )
  
  lam00 <- sum(hk*alp[-1])
  H1 <- as.matrix(hk)
  H2 <- diag( c((1-hk)/alp[-1]) )
  
  pi <- 0
  for (j in 0:1)
    for (k in 1:K) 
      pi <- pi + choose(K, k) * g.beta^k * (1-g.beta)^{K-k} * 
    (x.obs==j) * hk[j*K+k]
  
  V11 <- 1 - 1/alp[1]
  V31 <- as.matrix( -1/alp[1]*rep(1, 2*K) )
  V33 <- -1/alp[1]*rep(1, 2*K)%*%t(rep(1, 2*K)) - H2 + 
    1/nu*sum(1/pi^2)*H1%*%t(H1)
  
  V22 <- 0
  V23 <- 0
  V43 <- 0
  V24 <- 0
  V44 <- 0
  
  Uk <- NULL
  for(j in 0:1)
    for(k in 1:K){
      Uk <- cbind( Uk, choose(K, k) * g.beta^k * (1 - g.beta)^(K-k) *
                     (x.obs==j) - alp[j*K+k+1] )
    }
  
  for (i in 1:m) {
    pi.dot <- 0
    pi.ddot <- 0
    for(j in 0:1)
      for (k in 1:K){
        pi.dot <- pi.dot + choose(K, k) * 
          g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[j*K+k] *(x.obs[i]==j)*
          (k - K*g.beta[i]) * as.matrix(z.mat[i,]) 
        
        pi.ddot <- pi.ddot + choose(K, k) * 
          g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[j*K+k] *(x.obs[i]==j)*
          ( (k - K*g.beta[i])^2 - K*g.beta[i]*(1-g.beta[i]) ) * 
          as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) )
      }
    
    
    V22 <- V22 + ( pi.dot%*%t(pi.dot)/pi[i] - pi.ddot - 
                     K* g.beta[i] * (1-g.beta[i]) * pi[i] *
                     as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) ) )/pi[i]
    V23 <- V23 - pi.dot/pi[i]^2
    V43 <- V43 + Uk[i,]/pi[i]^2
    
    Uk.dot <- NULL
    for (j in 0:1)
      for(k in 1:K) Uk.dot <- cbind( Uk.dot, choose(K, k) * g.beta[i]^k * (1 - g.beta[i])^(K-k) * 
                                       (x.obs[i]==j) * (k - K*g.beta[i]) * z.mat[i, ] )
    
    V24 <- V24 - Uk.dot/pi[i] + pi.dot%*%t(as.matrix(Uk[i,]))/pi[i]^2
    V44 <- V44 + as.matrix(Uk[i,])%*%t(as.matrix(Uk[i,]))/pi[i]^2
  }
  V22 <- V22/nu
  V23 <- as.matrix(V23/nu)%*%t(H1)
  V43 <- lam00*( diag(1,2*K) - as.matrix(V43/nu)%*%t(H1) )
  V24 <- lam00*V24/nu
  V44 <- lam00^2*V44/nu
  
  inv.V44 <- solve(V44)
  # inv.V44 <- ginv(V44)
  
  W22 <- - V22 + V24%*%inv.V44%*%t(V24)
  W23 <- - V23 + V24%*%inv.V44%*%V43
  W33 <- - V33 + t(V43)%*%inv.V44%*%V43
  W <- rbind( cbind( -V11, matrix(0, 1, p), - t(V31) ),
              cbind( matrix(0, p, 1), W22, W23 ),
              cbind( -V31, t(W23), W33))
  
  ### variance estimates of MELE
  inv.W <- solve(W)
  # inv.W <- ginv(W)
  
  diag(inv.W)*c(nu, rep(1/nu,nrow(W) - 1))
}



#################################################################################
###  7. Calculate empirical log-likelihood under the hypothesis 
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
