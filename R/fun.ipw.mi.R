#################################################################################
###  Abundance estimation when covariates are subject to missing 
###        under discrete time capture-recapture models 
###  The inverse probability weighting and multiple imputation methods
#################################################################################

#################################################################################
###  1. Function to implement the inverse probability weighting method
#################################################################################

ipw.miss <- function(x.mis, z.mat, d.obs, d.mis, K){
  
  ##############################################################################
  ############   Preparation  ################
  
  p <- ncol(z.mat)
  y.obs <- z.mat[, p]
  d.all <- c( d.obs, d.mis) 
  n <- length(d.all)
  m <- length(d.obs)
  
  cov.all <- rbind(z.mat, cbind(rep(1, n-m), x.mis, rep(0, n-m)))
  ind.obs <- c( rep(1, m), rep(0, n-m) )
  
  if (p != 2) {
    
    x.obs <- z.mat[, p-1]
    x.all <- cov.all[,p-1]
    
  }
  
  ##############################################################################
  ############   Main functions  ################
  
  pi.yx.non.para.fun <- function(){
    
    out <- NULL
    for (i in 1:n) {
      
      if (p == 2) {
        
        n.obs <- sum( d.obs == d.all[i] )
        n.all <- sum( d.all == d.all[i] )
        
      }else{
        
        n.obs <- sum( d.obs == d.all[i] & x.obs == x.all[i] )
        n.all <- sum( d.all == d.all[i] & x.all == x.all[i] )
        
      }
      
      out <- c(out, n.obs/n.all )
      
    }
    
    out
  }
  
  pi.yx <- pi.yx.non.para.fun()
  
  pi.yx.non.na <- 0
  for (i in 1:n) pi.yx.non.na[i] <- ifelse( pi.yx[i] == 0, 1e-300, pi.yx[i] )
  
  like.beta.ipw <- function(beta){
    
    g.beta <- expit( z.mat %*% beta )
    phi.beta <- 1 - (1 - g.beta)^K
    - sum( (d.obs * log(g.beta + 1e-300) + (K - d.obs)*log(1 - g.beta + 1e-300) - log(phi.beta + 1e-300))/pi.yx[1:m] )
    
  }
  
  ###### IPW estimator of nu
  beta.est <- nlminb(rep(0,p), like.beta.ipw)$par
  gi <- plogis(cov.all%*%beta.est)
  p.star <- 1 - ( 1 - gi )^K
  nu.est <- sum( ind.obs/( pi.yx.non.na * p.star ) )
  
  
  ### Variance estimator of beta
  psi.i <- NULL
  psi.i.star <- matrix(0,n,p)
  kap.star <- rep(0, n)
  for (i in 1:n) psi.i <- rbind( psi.i, (d.all[i] - K*gi[i]/p.star[i]) * cov.all[i,] )
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
    cov.all[i,] %*% t(cov.all[i,]) / ( p.star[i]^2 * pi.yx[i] )
  
  var.beta <- solve(Gw) %*% ( t(Uwi) %*% Uwi ) %*% t(solve(Gw))
  
  ### Variance estimator of nu
  A.hat <- t(- t(cov.all) %*% ( ind.obs * K * gi * (1 - p.star) / (pi.yx.non.na * p.star^2) ) ) %*% t(solve(Gw))
  
  out <- ind.obs / (pi.yx.non.na * p.star) - kap.star * (ind.obs - pi.yx.non.na) / pi.yx.non.na +
    ( ind.obs * psi.i / pi.yx.non.na - (ind.obs - pi.yx.non.na) * psi.i.star) %*% t(A.hat)
  sigma2.nu <- sum(out^2) - nu.est
  
  list(nu.est = nu.est, sigma2.nu = sigma2.nu, 
       beta.est = beta.est, sigma2.beta =  diag(var.beta) )
}


#################################################################################
###  2. Function to implement the multiple imputation method
#################################################################################

im2.miss <- function(x.mis, z.mat, d.obs, d.mis, K, M){
  
  ##############################################################################
  ############   Preparation  ################
  
  p <- ncol(z.mat)
  d.all <- c( d.obs, d.mis) 
  n <- length(d.all)
  m <- length(d.obs)
  
  if (p == 2) {
    data <- cbind(d.all, 1, c(z.mat[, p], rep(0, n-m)), c(rep(1, m), rep(0, n-m)) )
  }else{
    x.all <-  c(z.mat[, p-1], x.mis)
    data <- cbind(d.all, 1, x.all, c(z.mat[, p], rep(0, n-m)), c(rep(1, m), rep(0, n-m)) )
    
  } 
  
  ##############################################################################
  ############   Main functions  ################
  
  MI.data.fun <- function(){
    
    mi.data <- matrix(0, n, M+1)  ###对Z填补后的10个数据集，第11列为缺失示性变量
    colnames(mi.data) <- c( paste0('Z',1:M), 'R')
    mi.data[,M+1] <- data[, p+2]
    for (i in 1:n) {
      
      if (data[i, p+2] == 1) {
        mi.data[i, 1:M] <- rep(data[i, p+1], M)
      }else{
        
        xx <- data[ d.all == d.all[i] & data[, p] == data[i, p] & data[, p+2] == 1, p+1]
        if (length(xx) != 0) {
          mi.data[i, 1:M] <- sample(xx, M, replace = T)
        }else{
          mi.data[i,1:M] <- rep( mean(data[1:m, p+1]), M)
        }
        
      }
    }
    
    mi.data
  }
  
  like.con.mi <- function(beta){
    
    like.i <- function(d, p.i, p.i.star) d*log(p.i + 1e-300) + (K - d)*log(1 - p.i + 1e-300) - log(p.i.star + 1e-300)
    
    out <- 0
    for (mm in 1:M) {
      
      if(p == 2){
        p.i <- expit(beta[1] + beta[2]*mi.data[,mm])
      }else{
        p.i <- expit(beta[1] + beta[2]*x.all + beta[3]*mi.data[,mm])
      }
      
      p.i.star <- 1 - (1 - p.i)^K
      out <- out + sum( like.i(d.all, p.i, p.i.star) )
    }
    - out/M
  }
  
  df.2 <- function(beta0){
    
    out <- 0
    for(mm in 1:M){
      
      if(p == 2){
        xx <- cbind(1, mi.data[,mm])
      }else{
        xx <- cbind(1, x.all, mi.data[,mm])
      }   
      
      p.i <- expit(xx%*%beta0)
      p.i.star <- 1 - (1-p.i)^K
      item <- - as.numeric( (K * p.i * (1-p.i) * p.i.star - K^2 * p.i^2 * (1-p.i.star) )/(p.i.star^2) )
      out <- out + t(xx) %*% diag(item) %*% xx
    }
    out/M
  }
  
  N.est.2 <- function(beta){
    
    p.i.star.mi <- 0
    
    for (mm in 1:M) {
      if(p == 2){
        p.i <- expit(beta[1] + beta[2]*mi.data[,mm])
      }else{
        p.i <- expit(beta[1] + beta[2]*x.all + beta[3]*mi.data[,mm])
      } 
      
      p.i.star <- 1 - (1 - p.i)^K
      p.i.star.mi <- p.i.star.mi + 1/p.i.star
    }
    sum(1/(M/p.i.star.mi))
  }
  
  var.est.2 <- function(beta){
    
    GG <- - df.2(beta)
    
    u.nu.i <- 0; u.ni <- 0; v.1 <- 0; v.3 <- 0; nu.est.mi <- NULL
    for (mm in 1:M) {
      
      if(p == 2){
        xx <- cbind(1, mi.data[,mm])
      }else{
        xx <- cbind(1, x.all, mi.data[,mm])
      }
      
      p.i <- expit(xx%*%beta)
      p.i.star <- 1 - (1 - p.i)^K
      u.nu.i <- u.nu.i + t(xx) %*% diag(as.numeric( (d.all - K*p.i/p.i.star)^2 )) %*% xx
      u.ni <- u.ni + t(xx) %*% (d.all - K*p.i/p.i.star) %*% t( t(xx) %*% (d.all - K*p.i/p.i.star) )
      v.1 <- v.1 + sum( (1-p.i.star)/p.i.star^2 ) 
      v.3 <- v.3 - t(xx) %*% ( K*(1 - p.i.star)*p.i/(p.i.star^2) )
      nu.est.mi <- c( nu.est.mi, sum(1/p.i.star) )
      
    }
    
    vv <- solve(GG) %*% (u.nu.i/M + (1+1/M) * u.ni / (M-1)) %*% solve(GG)
    
    c( diag(vv), v.1/M + (1+1/M) * sum( (nu.est.mi - est.mi[p+1])^2 )/(M-1) + t(v.3/M) %*% vv %*% (v.3/M) )
    
  }
  
  mi.data <- MI.data.fun()
  
  est.mi <- nlminb(rep(0, p), like.con.mi)$par
  est.mi <- c(est.mi, N.est.2(est.mi))
  est.var <- var.est.2( est.mi[1:p] )
  
  list(nu.est = est.mi[p+1], sigma2.nu = est.var[p+1], 
       beta.est = est.mi[1:p], sigma2.beta = est.var[1:p])
}
