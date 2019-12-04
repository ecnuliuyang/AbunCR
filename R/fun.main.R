#################################################################################
###   Main function for capture recapture model when covariates are subject to missing
###         Implement the empirical likelihood methods
#################################################################################

abun.opt <- function (d, K = NULL, x = NULL, y = NULL,
                      beta.initial = NULL,
                      CI = TRUE, level = 0.95,
                      SE = TRUE ) {

  if (is.null(K)) stop("please specify the number K")

   set.seed(321)
  if (is.null(beta.initial)) beta.initial <- runif( ncol(as.matrix(cbind(x,y))), -1, 1 )

  if ( is.null(y) ) {
    rt <- opt.el.cc ( d = d, K = K, x = x,
                      beta.initial = beta.initial, N0 = NULL )
  } else {

    rt <- opt.el.mar ( d = d, K = K, x = x, y = y,
                       beta.initial = beta.initial, N0 = NULL )

  }

  if ( CI ) {
    rt$n.big.ci <- abun.el.ci(rt, level = level)
  }
  rt$level <- level

  if ( SE ) {
    se <- abun.el.se(obj = rt)
    rt$n.big.se <- se$n.big.se
    rt$beta.se <- se$beta.se
  }

  class(rt) <- 'abun.opt'

  return(rt)

}


#################################################################################
###  Function to obtain the EL ratio confidence interval of N
#################################################################################

abun.el.ci <- function ( obj, level ) {

  ###  maximum empirical log-likelihood
  like.full <- obj$like
  n <- obj$n

  rn <- function (N0) {

    if (obj$class == 'abun.h')
      # (d, K, z, maxN, beta.initial, N0 )
      like.null <- opt.el.cc ( d = obj$d, K = obj$K, x = obj$x,
                               beta.initial = obj$beta.initial, N0 = N0)$like

    if (obj$class == 'abun.h.mar')
      like.null <-  opt.el.mar ( d = obj$d, K = obj$K, x = obj$x, y = obj$y,
                                 beta.initial = obj$beta.initial, N0 = N0 )$like

    2 * ( like.full - like.null ) - qchisq(level, 1)
  }


  hatN <- obj$n.big
  ntemp <- 5 * hatN
  while ( rn(ntemp) <= 0 )  ntemp <- ntemp*2

  if(ntemp > 1e5*n) stop('the likelihood is so flat that the upper limit is infinite')

  ci.upper <- uniroot( rn, c(hatN, ntemp), tol=0.01 )$root
  if ( rn(n) <= 0 ) {
    ci.lower <- n
  } else {
    ci.lower <- uniroot( rn, c(n, hatN), tol=0.01 )$root
  }

  ci <- c(ci.lower, ci.upper)

  return(ci)

}



#################################################################################
###  Function to calculate the Std. Error of estimates
#################################################################################
abun.el.se <- function (obj) {

  if (obj$class == 'abun.h')
    rt <- se.el.cc( obj )

  if (obj$class == 'abun.h.mar') {
    rt <- se.el.mar( obj )
  }

  return(rt)

}


#################################################################################
###  Function to show the results
#################################################################################

print.abun.opt <- function (obj) {

  cat ("\nCoefficients:\n")
  if ( is.null(obj$beta.se) ) {
    coefficients <- matrix( round(obj$beta, 2),1)
    rownames(coefficients) <- "Estimate"

  } else {
    beta <- matrix(round(obj$beta, 2),nrow=1)
    coefficients <- rbind( beta,
                           round(obj$beta.se[1:ncol(beta)], 2) )
    rownames(coefficients) <- c('Estimate','Std. Error')
  }
  colnames(coefficients) <- paste0('beta', 1:ncol(coefficients))

  print(coefficients)

  cat ("\nAbundance:\n")
  cat ("Estimate:", round(obj$n.big))
  if ( is.null(obj$n.big.se) ) {
    cat ("\nStd. Error: NULL")
  } else {
    cat ("\nStd. Error:", round(obj$n.big.se))
  }

  if ( is.null(obj$n.big.ci) ) {
    cat ("\n")
    cat (paste0(100*obj$level,"%"),
         "EL ratio CI: NULL \n")
  } else {
    cat ("\n")
    cat (paste0(100*obj$level,"%"),
         "EL ratio CI: [",
         paste(round(obj$n.big.ci), collapse=', '), "]\n")
  }

  cat("\nAIC of model:", obj$AIC,"\n\n")
}

