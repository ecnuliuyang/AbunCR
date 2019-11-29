#################################################################################
###   Main function for capture recapture model when covariates are subject to missing
###         Implement the empirical likelihood methods
#################################################################################

abun.opt <- function (formula, data, K = NULL, 
                      maxN = NULL, beta.initial = NULL,
                      CI = TRUE, level = 0.95, 
                      SE = TRUE ) {
  
  if (is.null(K)) stop("please specify the number K")
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action <- na.pass
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  xlev <- .getXlevels(mt, mf)
  
  d <- as.matrix( model.response(mf, "any") )
  x <- if (!is.empty.model(mt))
    model.matrix(mt, mf)
  
  if (is.null(beta.initial)) beta.initial <- runif( ncol(as.matrix(x)), -1, 1 )
  
  if ( all(!is.na(x)) ) {
    rt <- opt.el.cc ( d = d, K = K, z = as.matrix(x), maxN = NULL,
                      beta.initial = beta.initial, N0 = NULL )
  } else {
    
    if ( !is.matrix(x) ) {
      x.x = NULL
    } else {
      x.x <- x[, - unique(which(is.na(x), arr.ind = T)[,2])]
    }
    
    x.y <- x[, unique(which(is.na(x), arr.ind = T)[,2])]
    # print(head(x.x))
    # print(head(x.y))
    rt <- opt.el.mar ( d = d, K = K, x = x.x, y = x.y, 
                       maxN = maxN, beta.initial = beta.initial, N0 = NULL )
    
  }
  
  if ( CI ) {
    rt$n.big.ci <- abun.el.ci(rt, level = level)
  }
  rt$level <- level
  
  if ( SE )  rt$se <- abun.el.se(obj = rt)
  
  rt$call <- cl
  if (is.matrix(x))  rt$x.names <- colnames(x)
  else rt$x.names <- attr(mt,"term.labels")
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
      like.null <- opt.el.cc ( d = obj$d, K = obj$K, z = obj$z, maxN = obj$maxN,
                               beta.initial = obj$beta.initial, N0 = N0)$like
    
    if (obj$class == 'abun.h.mar')
      like.null <-  opt.el.mar ( d = obj$d, K = obj$K, x = obj$x, y = obj$y, maxN = obj$maxN, 
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
  
  cat("\nCall:\n", paste(deparse(obj$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")   
  
  cat ("Coefficients:\n")
  if ( is.null(obj$se) ) {
    coefficients <- matrix( round(obj$beta, 2),1)  
    rownames(coefficients) <- "Estimate"
    
  } else {
    coefficients <- rbind( round(obj$beta, 2),
                           round(obj$se$beta.se, 2) )
    rownames(coefficients) <- c('Estimate','Std. Error')
  }
  colnames(coefficients) <- paste0('beta', 1:ncol(coefficients))
  
  print(coefficients)
  
  cat ("\nAbundance:\n")
  cat ("Estimate:", round(obj$n.big))
  if ( is.null(obj$se) ) {
    cat ("\nStd. Error: NULL")
  } else {
    cat ("\nStd. Error:", round(obj$se$n.big.se))
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


