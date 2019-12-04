#################################################################################
### Empirical likelihood method in the presence of missing data
#################################################################################

#################################################################################
###  1. Calculate Ak for a given beta
#################################################################################

ak.fun <- function (g.beta, x.obs, lev.mis, K, p) {
  Ak <- (1 - g.beta)^K
  for (i in 1:nrow(lev.mis))
    Ak <- cbind( Ak, choose(K, lev.mis[i,1]) *
                   g.beta^lev.mis[i,1] * (1 - g.beta)^(K-lev.mis[i,1]) *
                   apply( x.obs, 1, function(x) all(x == lev.mis[i,-1]) ) )
  Ak
}


#################################################################################
###  2. Minimize l_3(N, beta, lambda) with respect to lambda
#################################################################################

like3.mar <- function (g.beta, x.obs, lev.mis, n, m, mk, K, n.big, p) {

  logg <- function (t, eps = 1e-5) {
    ts <- t*(t>eps)+1*(t<=eps)
    log(ts) * (t>eps) + (log(eps) - 1.5 + 2*t/eps - (t/eps)^2*0.5) * (t<=eps)
  }

  el.lam <- function (lam)
    - sum( logg( Ak%*%lam + n.big/m) ) - sum( c(n.big - n, mk) * log(-lam) )

  lambda0 <- c( - (n.big - n)/m, - (mk/m) )
  Ak <- ak.fun(g.beta, x.obs, lev.mis, K, p)

  nlminb(lambda0 - 1e-5, el.lam, upper = lambda0 - 1e-20)

}


#################################################################################
###  3. Maximize l_1(N) + min_lambda l_3(N, beta, lambda) with respect to N
###  for the given beta
#################################################################################

like13.mar <- function (g.beta, x.obs, lev.mis, n, m, mk, K, p, N0 = NULL) {

  like.n.big <- function(n.big){

    like1 <- sum( log( (n.big - n + 1):n.big ) ) + (n.big - n) * log( (n.big - n)/m + 1e-300 )
    like3 <- like3.mar(g.beta, x.obs, lev.mis, n, m, mk, K, n.big, p)$objective
    - (like1 + like3)

  }

  if (is.null(N0)) {
    out <- optimize(like.n.big, interval=c(n, 100*n), tol=0.01)
  } else {
    out <- like.n.big(N0)
  }
  out

}


#################################################################################
###  4. Implement empirical likelihood method in the presence of missing data
#################################################################################

opt.el.mar <- function ( d, K, x, y,
                          beta.initial, N0 = NULL) {
  d <- as.numeric(d)
  n <- length(d)

  if( is.null(x) ) x.mat <- matrix(1, nrow=n)
  else  x.mat <- as.matrix(x)

  y.mat <- as.matrix(y)
  na.loc <- apply( y.mat, 1, function(y) any(is.na(y)) )
  m <- sum(!na.loc)
  d.obs <- d[!na.loc]
  d.mis <- d[na.loc]
  x.mis <- as.matrix(x.mat[na.loc,])
  x.obs <- as.matrix(x.mat[!na.loc,])
  y.obs <- as.matrix(y.mat[!na.loc,])

  lev.mis <- as.data.frame( cbind(d.mis, x.mis) )
  group.sum <- aggregate(rep(1,n-m),
                         by = as.list(lev.mis), sum)

  lev.mis <- group.sum[, 1:(ncol(group.sum) - 1)]
  colnames(lev.mis) <- c('d', paste0('x',1:ncol(x.obs)))
  lev.mis <- as.matrix(lev.mis)
  mk <- group.sum[, ncol(group.sum)]


  ### z.obs = z.mat
  z.obs = cbind(x.obs, y.obs)
  p <- ncol(z.obs)
  z.mat <- as.matrix(z.obs[, beta.initial!=0])

  el.beta <- function(beta){
    beta <- as.matrix(beta)
    g.beta <- as.numeric( plogis( z.mat %*% beta ) )
    like2.mar <- sum( d.obs * log(g.beta + 1e-300) + (K - d.obs) * log( 1 - g.beta + 1e-300 ) )
    temp <- like13.mar(g.beta, x.obs, lev.mis, n, m, mk, K, p, N0 = N0)

    if (is.null(N0)) {
      out <- - ( like2.mar - temp$objective )
    } else {
      out <- - ( like2.mar - temp )
    }
    out
  }


  out <- nlminb(beta.initial[beta.initial!=0], el.beta)
  like <- - out$objective
  beta.est <- as.matrix(out$par)
  g.beta <- as.numeric( plogis(z.mat%*%  beta.est ) )

  ##############  EL estimator of n.big  ##############
  out13 <- like13.mar(g.beta, x.obs, lev.mis, n, m, mk, K, p, N0 = N0)
  if (is.null(N0)) n.big.est <- out13$minimum
  else n.big.est <- N0

  ############## EL estimator of lam, alp and pi's ############
  out3 <- like3.mar(g.beta, x.obs, lev.mis, n, m, mk, K, n.big.est, p)
  lam.est <- out3$par
  alp.est <- - c(n.big.est - n, mk)/(m*lam.est)
  Ak <- ak.fun(g.beta, x.obs, lev.mis, K, p)
  prob.est <- as.numeric( 1/(m * (1 + ( Ak - matrix(alp.est, m, length(alp.est), byrow = T) )
                                  %*%lam.est)) )

  ###  Maximum emiprical log-likelihood & AIC  ##############
  log.like <- sum(log((n.big.est - n+1):n.big.est)) - sum(log(1:n)) +
    sum( c(n.big.est - n, mk) * log(alp.est + 1e-300) ) +
    sum( log(prob.est + 1e-300) ) +
    sum( d.obs * log(g.beta + 1e-300) +
           (K - d.obs) * log( 1 - g.beta + 1e-300 ) )
  case1 <- length(alp.est) ==
    K*prod(apply(x.mat, 2, function(x)
      length(unique(x)))) + 1

  if (case1) num.pars <- length(beta.est) + length(alp.est)
  else num.pars <- length(beta.est) + length(alp.est) + 1

  AIC <- 2*(-log.like + num.pars)

  rt <- list( n.big = n.big.est,
              beta = beta.est,
              alpha = alp.est,
              like = log.like, AIC = AIC,
              prob = prob.est,
              class = 'abun.h.mar',
              d = d, K = K, x = x, y = y,
              beta.initial=beta.initial,
              n=n, m=m)
  return(rt)

}


#################################################################################
###  5. Calculate the Std. Error of estimates
#################################################################################

se.el.mar <- function(obj) {

  d <- as.numeric(obj$d)
  K <- obj$K
  x <- obj$x
  y.mat <- as.matrix(obj$y)
  beta.initial <- obj$beta.initial
  beta <- as.matrix( obj$beta )
  p <- length(beta)
  n.big <- obj$n.big

  na.loc <- apply(y.mat, 1, function(y) any(is.na(y)))
  n <- length(na.loc)
  m <- sum(!na.loc)

  d.obs <- d[!na.loc]
  d.mis <- d[na.loc]
  y.obs <- as.matrix(y.mat[!na.loc,])

  if (is.null(x)) x.mat <- matrix(1, n, 1)
  else x.mat <- as.matrix(x)

  x.mis <- as.matrix(x.mat[na.loc,])
  x.obs <- as.matrix(x.mat[!na.loc,])

  lev.mis <- as.data.frame( cbind(d.mis, x.mis) )
  group.sum <- aggregate(rep(1,n-m),
                         by = as.list(lev.mis), sum)
  lev.mis <- group.sum[, 1:(ncol(group.sum) - 1)]
  colnames(lev.mis) <- c('d', paste0('x',1:ncol(x.obs)))
  lev.mis <- as.matrix(lev.mis)

  mk <- group.sum[, ncol(group.sum)]


  ### z.obs = z.mat
  z.obs = as.matrix( cbind(x.obs, y.obs) )
  z.mat <- as.matrix( z.obs[, beta.initial!=0] )
  g.beta <-as.numeric( plogis(z.mat%*%beta) )
  alp0 <- obj$alpha[1]
  pi <- obj$prob

  if (is.null(x)) {
    x.mat <- matrix(1, n, 1)
    hk = temp.alpha <- NULL
    for (k in 1:K){
      hk <- c( hk, sum( d.obs == k )/sum( c(d.obs, d.mis) == k ) )
      temp.alpha <- c(temp.alpha, choose(K, k) * sum( g.beta^k * (1 - g.beta)^{K - k} * pi))
    }

  } else {
    x.mat <- as.matrix(x)
    x.unique <- as.matrix( unique(x.mat) )
    hk = temp.alpha <- matrix(NA, nrow(x.unique), K)
    for (i in 1:nrow(x.unique)) {
      for (k in 1:K){
        hk[i,k] <- sum( d.obs == k & apply( x.obs, 1, function(x) all(x == x.unique[i,]) ) )/
                   sum( c(d.obs, d.mis) == k&apply( rbind(x.obs,x.mis), 1, function(x) all(x == x.unique[i,]) ) )
        temp.alpha[i,k] <- choose(K, k) *
          sum( g.beta^k * (1 - g.beta)^{K - k} *
                 apply( x.obs, 1, function(x) all(x ==  x.unique[i,]) ) * pi)
      }
    }

  }

  hk[is.na(hk)] <- 1
  Ak <- ak.fun(g.beta, x.obs, lev.mis, K, p)

  lam00 <- sum(hk*temp.alpha)

  hk.pos <- NULL
  alp.pos <- NULL
  for (ii in 1:nrow(lev.mis))
    for (j in 1:nrow(x.unique))
      for (k in 1:K) {
        if( lev.mis[ii,1] == k & all(lev.mis[ii,-1] == x.unique[j,]) ){
          hk.pos <- c(hk.pos, hk[j,k])
          alp.pos <- c(alp.pos, temp.alpha[j,k])

        }
      }

  Uk <- (Ak - matrix(c(alp0, alp.pos), nrow(Ak), ncol(Ak), byrow = TRUE))

  pis <- 0
  for (ii in 1:nrow(x.unique))
    for (k in 1:K)
      pis <- pis + choose(K, k) * g.beta^k * (1-g.beta)^{K-k} *
    apply( x.obs, 1, function(x) all(x == x.unique[ii,])) * hk[ii, k]


  V11 <- 1 - 1/alp0

  case1 <- length(obj$alpha) ==
    K*prod(apply(x.mat, 2, function(x)
      length(unique(x)))) + 1
  if ( case1 ) {

    H1 <- as.matrix( hk.pos )
    H2 <- diag( (1-hk.pos)/alp.pos )
    V31 <- diag( - 1/alp0, length(mk), length(mk) )
    V33 <- - 1/alp0 * matrix(1, nrow(H2s), ncol(H2s) ) - H2s +
      1/n.big*sum(1/pis^2)*H1s%*%t(H1s)

    Uk <- Uk[,-1]
  } else {

    H1s <- as.matrix( c(1, 1-hk.pos) )
    H2s <- diag( c(1/alp0, (1-hk.pos)/alp.pos) )
    V31 <- as.matrix( c(1/alp0, rep(0, length(mk))) )
    V33 <- - H2s + 1/n.big*sum(1/pis^2)*H1s%*%t(H1s)

  }


  V22 <- 0
  V23 <- 0
  V43 <- 0
  V24 <- 0
  V44 <- 0


  for (i in 1:m) {
    pi.dot <- 0
    pi.ddot <- 0
    for (j in 1:nrow(x.unique))
      for (k in 1:K) {
        pi.dot <- pi.dot + choose(K, k) *
          g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[j,k] *
          all(x.obs[i,] == x.unique[j,]) *
          (k - K*g.beta[i]) * as.matrix(z.mat[i,])

        pi.ddot <- pi.ddot + choose(K, k) *
          g.beta[i]^k * (1-g.beta[i])^{K-k} * hk[j,k] *
          all(x.obs[i,] == x.unique[j,]) *
          ( (k - K*g.beta[i])^2 - K*g.beta[i]*(1-g.beta[i]) ) *
          as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) )
      }



    V22 <- V22 + ( pi.dot%*%t(pi.dot)/pis[i] - pi.ddot -
                     K* g.beta[i] * (1-g.beta[i]) * pis[i] *
                     as.matrix(z.mat[i,])%*%t( as.matrix(z.mat[i,]) ) )/pis[i]
    if (case1)  {
      V23 <- V23 - pi.dot/pis[i]^2
      V43 <- V43 - Uk[i,]/pis[i]^2

    } else {
      V23 <- V23 + pi.dot/pis[i]^2
      V43 <- V43 + Uk[i,]/pis[i]^2

    }

    if (case1) Uk.dot <- NULL
    else Uk.dot <- (1 - g.beta[i])^K *
      (0 - K*g.beta[i]) * z.mat[i, ]


    for (ii in 1:nrow(lev.mis))
      for (j in 1:nrow(x.unique))
        for (k in 1:K) {
          if ( lev.mis[ii,1] == k & all(lev.mis[ii,-1] == x.unique[j,]) )
             Uk.dot <- cbind( Uk.dot, choose(K, k) * g.beta[i]^k * (1-g.beta[i])^(K-k)*
                               all(x.obs[i,] == x.unique[j,]) *
                               (k - K*g.beta[i]) * z.mat[i, ] )
          }

    V24 <- V24 - Uk.dot/pis[i] + pi.dot%*%t(as.matrix(Uk[i,]))/pis[i]^2
    V44 <- V44 + as.matrix(Uk[i,])%*%t(as.matrix(Uk[i,]))/pis[i]^2
  }

  V22 <- V22/n.big
  V23 <- as.matrix(V23/n.big)%*%t(H1s)
  V43 <- lam00*( diag(1,length(H1s)) + as.matrix(V43/n.big)%*%t(H1s) )
  V24 <- lam00*V24/n.big
  V44 <- lam00^2*V44/n.big
  inv.V44 <- solve(V44 - diag(1e-300, nrow(V44), ncol(V44)))

  W22 <- - V22 + V24%*%inv.V44%*%t(V24)
  W23 <- - V23 + V24%*%inv.V44%*%V43
  W33 <- - V33 + t(V43)%*%inv.V44%*%V43
  Ws <- rbind( cbind( -V11, matrix(0, 1, p), - t(V31) ),
               cbind( matrix(0, p, 1), W22, W23 ),
               cbind( -V31, t(W23), W33) )

  ### variance estimates of MELE
  inv.Ws <- solve(Ws - diag(1e-300, nrow(Ws), ncol(Ws)))
  # se <- sqrt( diag(inv.Ws)[1:(length(beta)+1)]*c(n.big, rep(1/n.big, length(beta))) )
  se <- sqrt( diag(inv.Ws) * c(n.big, rep(1/n.big, ncol(inv.Ws)-1)) )

  rt <- list()
  rt$n.big.se <- se[1]
  rt$beta.se <- se[2:(1+ length(beta))]
  return(rt)
}


