#################################################################################
###  Abundance estimation when covariates are subject to missing 
###        under discrete time capture-recapture models 
###  The naive empirical likelihood method
#################################################################################

#################################################################################
###  1. Minimize   - sum.i  log(1+lambda t.i)  with respect to lambda                      
#################################################################################

el.lagrange <- function(t){
  
  emplik <- function(lam)  - sum( logg( 1 + lam*t ) )
  
  optimize(f=emplik, interval=c(-1000, 1000))
}


#################################################################################
###  2. Minimize  - h.1(N,alpha) with respect to N for the given alpha    
#################################################################################

slike.part1 <- function(m, alpha){ 
  
  fgamma <- function(nu) - sum( log((nu - m + 1):nu) ) - (nu - m) * log(alpha)  
  
  optimize(fgamma, interval=c(m, 100*m), tol=0.01)
}

#################################################################################
###  3. Minimize -h.23(beta,alpha) with respect to beta for the given alpha 
#################################################################################

slike.part23 <- function(z.mat, d.obs, K, beta.init, alpha){
  
  fun23 <- function(beta){
    
    g.beta <- expit(z.mat%*%beta)
    phi <- (1 - g.beta)^K
    result.lagrange <- el.lagrange( phi - alpha )
    
    tt <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log(1 - g.beta + 1e-300) )
    
    - ( result.lagrange$objective + tt )
    
  }
  
  gr <- function(beta){
    
    g.beta <- expit(z.mat%*%beta)
    phi <- (1-g.beta)^K
    result.lagrange <- el.lagrange( phi - alpha )
    lam <- result.lagrange$minimum
    temp <- ( lam * alpha - 1 ) / ( 1 + lam * (phi - alpha) ) * K * g.beta
    
    - t(z.mat) %*% (temp + d.obs)
    
  }
  
  nlminb(beta.init, fun23, gradient=gr)
  
}


#################################################################################
###  4. Calculate empirical likelihood function for the given N.0  
#################################################################################
sel.null <- function(z.mat, d.obs, K, beta.init, nu0){
  
  m <- length(d.obs)
  falpha.null <- function(alpha){
    
    tmp1 <- - sum( log((nu0 - m + 1):nu0) ) - (nu0 - m) * log(alpha)
    tmp2 <- slike.part23(z.mat, d.obs, K, beta.init, alpha)$objective
    tmp1 + tmp2
  }
  
  eps <- 1e-5
  out <- optimize(falpha.null, interval=c(eps, 1-eps), tol=0.0001)
  - out$objective
}


#################################################################################
###  5. Implement empirical likelihood method    
#################################################################################
el.cc <- function(z.mat, d.obs, K, beta.init, level){
  
  m <- length(d.obs)
  falpha <- function(alpha){
    
    tmp1 <- slike.part1(m, alpha)$objective
    tmp2 <- slike.part23(z.mat, d.obs, K, beta.init, alpha)$objective
    tmp1 + tmp2
  }
  
  eps <- 1e-6
  out <- optimize(falpha, interval=c(eps, 1-eps), tol=0.0001)
  like <- - out$objective
  alpha.est <- out$minimum
  
  ##############  EL estimator of N  ##############
  
  out1 <- slike.part1(m, alpha.est)
  nu.est <- out1$minimum
  
  #################################################
  ############## EL estimator of beta  ############
  
  out23 <- slike.part23(z.mat, d.obs, K, beta.init, alpha.est)
  beta.est <- out23$par
  
  #################################################
  
  
  ###########   probability weights of EL  ##############
  
  g.beta <- expit(z.mat %*% beta.est)
  phi <- (1 - g.beta)^K
  lam <- el.lagrange(phi - alpha.est)$minimum
  prob.est <- 1 / ( m*(1+lam*(phi-alpha.est)) )
  
  ###########      Test of convergence     ##############
  converge <- ifelse( abs(1 - sum(prob.est) ) < 1e-2 & all(prob.est > 0) & all(prob.est < 1), 0, 1 )
  
  ###########   ELR based CI  ##############
  
  rn <- function(nu){
    
    like.full <- like
    like.null <- sel.null(z.mat, d.obs, K, beta.est, nu)
    2 * ( like.full - like.null) - qchisq(level, 1)
    
  }
  
  ntemp <- nu.est
  while(rn(ntemp) <= 0) ntemp <- ntemp*2
  
  ci.upper <- uniroot(rn, c(nu.est, ntemp))$root
  if (rn(m) <= 0) {
    ci.lower <- m
  }else{
    ci.lower <- uniroot(rn, c(m, nu.est))$root
  }
  
  ci.est <- c(ci.lower, ci.upper)

  
  list(nu.est = nu.est, ci.est = ci.est, alpha.est = alpha.est, beta.est = beta.est, converge = converge)
  
}


