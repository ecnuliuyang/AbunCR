#################################################################################
###  Abundance estimation when covariates are subject to missing 
###        under discrete time capture-recapture models 
###  The naive empirical likelihood method
#################################################################################

#################################################################################
###  1. Minimize   - sum.i  log(1+lambda t.i)  with respect to lambda                      
#################################################################################

el.lagrange <- function(t){
  
  ### Function to approximate log(t)
  logg <- function(t){
    
    eps <- 1e-5
    ts <- t*(t>eps)+1*(t<=eps)
    log(ts) * (t>eps) + (log(eps) - 1.5 + 2*t/eps - (t/eps)^2*0.5) * (t<=eps)
    
  }
  
  emplik <- function(lam)  - sum( logg( 1 + lam*t ) )
  
  optimize(f=emplik, interval=c(-1e3, 1e3))
}


#################################################################################
###  2. Minimize  - h.1(N,alpha) with respect to N for the given alpha    
#################################################################################

slike.part1 <- function(n, alpha){
  
  fgamma <- function(nt) - sum( log((nt - n + 1):nt) ) - (nt - n) * log(alpha+1e-20)  
  
  optimize(fgamma, interval=c(n, 100*n), tol=0.01)
  
}

#################################################################################
###  3. Minimize -h.23(beta,alpha) with respect to beta for the given alpha 
#################################################################################

slike.part23 <- function(z.mat, d.obs, K, beta.init, alpha){
  
  z.mat <- as.matrix(z.mat)
  
  fun23 <- function(beta){
    
    g.beta <- as.numeric( plogis(z.mat%*%beta) )
    phi <- (1 - g.beta)^K
    result.lagrange <- el.lagrange( phi - alpha )
    
    tt <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log(1 - g.beta + 1e-300) )
    
    - ( result.lagrange$objective + tt )
    
  }
  
  gr <- function(beta){
    
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
###  4. Calculate empirical likelihood function for the given N.0  
#################################################################################
sel.null <- function(z.mat, d.obs, K, beta.init, n.big0){
  
  n <- length(d.obs)
  falpha.null <- function(alpha){
    
    tmp1 <- - sum( log((n.big0 - n + 1):n.big0) ) - (n.big0 - n) * log(alpha)
    tmp2 <- slike.part23(z.mat, d.obs, K, beta.init, alpha)$objective
    tmp1 + tmp2
  }
  
  eps <- 1e-5
  out <- optimize(falpha.null, interval=c(eps, 1-eps), tol=eps)
  - out$objective
}


#################################################################################
###  5. Implement empirical likelihood method with no missing data   
#################################################################################
###  note: first column of z.mat is a vector with elements being 1.
el.opt.cc <- function(z, d, K, CI = TRUE, level = .95,
                   beta.initial = rep(0, ncol(as.matrix(z)))){
  
  z.mat <- as.matrix(z)
  d.obs <- as.numeric(d)
  n <- length(d.obs)
  falpha <- function(alpha){
    
    tmp1 <- slike.part1(n, alpha)$objective
    tmp2 <- slike.part23(z.mat, d.obs, K, beta.initial, alpha)$objective
    tmp1 + tmp2
  }
  
  eps <- 1e-5
  out <- optimize(falpha, interval=c(eps, 1-eps), tol=eps)
  like <- - out$objective
  alpha.est <- out$minimum
  
  ##############  EL estimator of N  ##############
  
  out1 <- slike.part1(n, alpha.est)
  n.big.est <- out1$minimum
  
  #################################################
  ############## EL estimator of beta  ############
  
  out23 <- slike.part23(z.mat, d.obs, K, beta.initial, alpha.est)
  beta.est <- out23$par
  
  #################################################
  
  
  ###########   probability weights of EL  ##############
  
  g.beta <- as.numeric( plogis(z.mat %*% beta.est) )
  phi <- (1 - g.beta)^K
  lam <- el.lagrange(phi - alpha.est)$minimum
  prob.est <- 1 / ( n*(1+lam*(phi-alpha.est)) )
  
  ###########      Test of convergence     ##############
  converge <- ifelse( abs(1 - sum(prob.est) ) < 1e-2 & all(prob.est > 0) & all(prob.est < 1), 0, 1 )
  
  #### Maximum empirical log-likelihood and AIC #######
  log.like <- sum( log((n.big.est - n + 1):n.big.est) ) - sum(log(1:n)) + 
    (n.big.est - n) * log(alpha.est) + sum(log(prob.est)) +
    sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log(1 - g.beta + 1e-300) )
  AIC <- 2*(-log.like + length(beta.est) + 1 + 1)
  
  ###########   ELR based CI  ##############
  
  rn <- function(nu){
    
    like.full <- like
    like.null <- sel.null(z.mat, d.obs, K, beta.est, nu)
    2 * ( like.full - like.null) - qchisq(level, 1)
    
  }
  
  if (CI) {
    ntemp <- n.big.est
    while(rn(ntemp) <= 0) ntemp <- ntemp*2
    
    ci.upper <- uniroot(rn, c(n.big.est, ntemp))$root
    if (rn(n) <= 0) {
      ci.lower <- n
    }else{
      ci.lower <- uniroot(rn, c(n, n.big.est))$root
    }
    
    ci.est <- c(ci.lower, ci.upper)
  } else ci.est <- NULL
  
  
  ###########   Variance estimates  ##############
  dim_q=ncol(z.mat)
  V23 = rep(0, dim_q)
  V22= matrix(0,dim_q,dim_q)
  
  beta=beta.est
  phi=(1-g.beta)^K
  n.big.est = n.big.est
  phistar=sum(1/(1-phi+1e-300)^2)/n.big.est
  
  V23.coef=(phi/(1-phi+1e-300)^2)*K*g.beta
  V22.coef=(d.obs-K*g.beta/(1-phi+1e-300))^2
  V22.part2=z.mat
  
  for(i in 1:n)
  {
    
    V23=V23+V23.coef[i]*z.mat[i,]
    qxi=matrix(z.mat[i,],ncol=1)
    
    V22=V22+V22.coef[i]*qxi%*%t(qxi)
    
  }
  
  V23=V23/n.big.est;
  V22=-V22/n.big.est
  V23=matrix(V23,ncol=1)
  V32=t(V23)
  V22.inv=solve(V22-diag(1e-300,dim_q,dim_q));
  tmp = V32%*%V22.inv%*%V23
  tmp=as.numeric(tmp)
  
  sigma_hat =  phistar - 1 - tmp
  
  list(n.big.est = n.big.est, ci.est = ci.est, 
       alpha.est = alpha.est, beta.est = beta.est,
       like = log.like, AIC = AIC,
       converge = converge,
       se.nu = sqrt(sigma_hat*n.big.est), 
       se.beta = sqrt(diag(- V22.inv)/n.big.est))
  
}


