#################################################################################
### Empirical likelihood method in the absence of missing data
#################################################################################

#################################################################################
###  1. Minimize   - sum.i  log(1+lambda t.i)  with respect to lambda                      
#################################################################################

el.lagrange <- function(t){
  
  ### Function to approximate log(t)
  logg <- function (t, eps=1e-5) {
    ts <- t*(t>eps)+1*(t<=eps)
    log(ts) * (t>eps) + (log(eps) - 1.5 + 2*t/eps - (t/eps)^2*0.5) * (t<=eps)
  }
  
  emplik <- function(lam)  - sum( logg( 1 + lam*t ) )
  
  optimize(f=emplik, interval=c(-1e3, 1e3))
  
}


#################################################################################
###  2. Minimize  - h.1(N,alpha) with respect to N for the given alpha    
#################################################################################

slike.part1 <- function(n, alpha, N0){
  
  f.nt <- function(nt) - sum( log((nt - n + 1):nt) ) - (nt - n) * log(alpha+1e-20)  
  
  if (is.null(N0)) 
    out <- optimize(f.nt, interval=c(n, 100*n), tol=0.01)
  else out <- f.nt(N0)
  out
}

#################################################################################
###  3. Minimize -h.23(beta,alpha) with respect to beta for the given alpha 
#################################################################################

slike.part23 <- function(z.mat, d.obs, K, beta.init, alpha){
  
  fun23 <- function(beta){
    
    beta <- as.matrix(beta)
    g.beta <- as.numeric( plogis(z.mat%*%beta) )
    phi <- (1 - g.beta)^K
    result.lagrange <- el.lagrange( phi - alpha )
    
    tt <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log(1 - g.beta + 1e-300) )
    
    - ( result.lagrange$objective + tt )
    
  }
  
  gr <- function(beta){
    beta <- as.matrix(beta)
    g.beta <- as.numeric( plogis(z.mat%*%beta) )
    phi <- (1-g.beta)^K
    result.lagrange <- el.lagrange( phi - alpha )
    lam <- result.lagrange$minimum
    temp <- ( lam * alpha - 1 ) / ( 1 + lam * (phi - alpha) ) * K * g.beta
    
    - t(z.mat) %*% (temp + d.obs)
    
  }
  
  nlminb(beta.init, fun23, gradient=gr)
  
}


#################################################################################
###  4. Implement empirical likelihood method with no missing data   
#################################################################################
###  note: first column of z.mat is all one.

opt.el.cc <- function (d, K, x, beta.initial, N0 ) {

  z.mat <- as.matrix( as.matrix(x)[, beta.initial!=0] )
  d.obs <- as.numeric(d)
  n <- length(d.obs)
  falpha <- function(alpha) {
    if (is.null(N0)) 
      tmp1 <- slike.part1(n, alpha, N0)$objective
    else
      tmp1 <- slike.part1(n, alpha, N0)
      
    tmp2 <- slike.part23(z.mat, d.obs, K, beta.initial[beta.initial!=0], alpha)$objective
    tmp1 + tmp2
  }
  
  eps <- 1e-5
  out <- optimize(falpha, interval=c(eps, 1-eps), tol=eps)
  # like <- - out$objective
  alpha.est <- out$minimum
  
  ##############  EL estimator of N  ##############
  if (is.null(N0)) {
    out1 <- slike.part1(n, alpha.est, N0)
    n.big.est <- out1$minimum
  } else n.big.est <- N0
 
  
  #################################################
  ############## EL estimator of beta  ############
  
  out23 <- slike.part23(z.mat, d.obs, K, beta.initial[beta.initial!=0], alpha.est)
  beta.est <- as.matrix(out23$par)
  
  #################################################
  ###########   probability weights of EL  ########
  
  g.beta <- as.numeric( plogis(z.mat %*% beta.est) )
  phi <- (1 - g.beta)^K
  lam <- el.lagrange(phi - alpha.est)$minimum
  prob.est <- 1 / ( n*(1+lam*(phi-alpha.est)) )
  
  #### Maximum empirical log-likelihood and AIC #######
  log.like <- sum( log((n.big.est - n + 1):n.big.est) ) - sum(log(1:n)) + 
    (n.big.est - n) * log(alpha.est) + sum(log(prob.est)) +
    sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log(1 - g.beta + 1e-300) )
  AIC <- 2*(-log.like + length(beta.est) + 2)

  rt <- list ( n.big = n.big.est, 
               beta = beta.est,
               alpha = alpha.est, 
               like = log.like, AIC = AIC,
               prob = prob.est, class = 'abun.h',
               d = d, K = K, x = x, 
               beta.initial = beta.initial,
               n = n)
  return(rt)
  
}

#################################################################################
###  5. Calculate the Stand. Error
#################################################################################

se.el.cc <- function(obj) {
  
  d.obs <- as.numeric(obj$d)
  n <- length(d.obs)
  z.mat <- as.matrix( as.matrix(obj$x)[,obj$beta.initial!=0] )

  dim_q <- ncol(z.mat)
  V23 <- rep(0, dim_q)
  V22 <- matrix(0, dim_q, dim_q)
  
  n.big <- obj$n.big
  K <- obj$K
  beta <- as.matrix(obj$beta)
  g.beta <- as.numeric( plogis(z.mat %*% beta) )
  phi <- (1-g.beta)^K
  phistar <- sum(1/(1-phi+1e-300)^2)/n.big
  
  V23.coef <- (phi/(1-phi+1e-300)^2)*K*g.beta
  V22.coef <- (d.obs-K*g.beta/(1-phi+1e-300))^2
  V22.part2 <- z.mat
  
  for(i in 1:n) {
    
    V23 <- V23+V23.coef[i]*z.mat[i,]
    qxi <- as.matrix(z.mat[i,])
    
    V22 <- V22+V22.coef[i]*qxi%*%t(qxi)
    
  }
  
  V23 <- V23/n.big
  V22 <- -V22/n.big
  V23 <- matrix(V23,ncol=1)
  V32 <- t(V23)
  V22.inv <- solve(V22-diag(1e-300,dim_q,dim_q))
  tmp <- V32%*%V22.inv%*%V23
  tmp <- as.numeric(tmp)
  
  sigma_hat  <-   phistar - 1 - tmp
  
  rt <- list()
  rt$n.big.se <- sqrt( sigma_hat*n.big )
  rt$beta.se <- sqrt( diag(- V22.inv)/n.big )
  return(rt)
}



