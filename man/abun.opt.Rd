
\name{abun.opt}
\alias{abun.opt}


\title{Maximum empirical likelihood (EL) abundance estimation under the Huggins-Alho (H-A) model}


\description{
This function is used to calculate the maximum EL estimate and EL ratio confidence interval of the abundance under the H-A model, as well as the EL estimates of parameters in the H-A model. Thier standard errors are also provided for reference.
}


\usage{

abun.opt(d, K = NULL, x = NULL, y = NULL,
  CI = TRUE, level = 0.95, SE = TRUE,
  beta.initial = NULL)

}



\arguments{

\item{d}{a vector, the number of times of being captured.}

\item{K}{a number, the number of capture occasions.}

\item{x}{a vector, matrix or data.frame, stands for individual covariates without missingness.}

\item{y}{a vector, matrix or data.frame, stands for individual covariates  with missing values (\code{NA}).}

\item{CI}{logistic, indiates whether the EL ratio confidence interval of the abundance is given (TRUE, default) or not (FALSE).}

\item{level}{a number, the nominal confidence level (default: 0.95).}

\item{SE}{logistic, indiates whether the Std. Error is given (TRUE, default) or not (FALSE).}

\item{beta.initial}{a vector, with the same length as the parameters in H-A model. 0 means that the coefficient is fixed to zero, and non-zero values are seen as initial values when optimizing.
}

}

\value{

A list containing the following main components:

\item{n.big}{the maximum EL estimate of the abundance.}

\item{n.big.se}{Std. Error of \code{n.big}.}

\item{n.big.ci}{the EL ratio confidence interval of the abundance.}

\item{beta}{the maximum EL estimates of the parameters in the H-A model.}

\item{beta.se}{Std. Error of \code{beta}.}

\item{like}{the maximum EL value.}

\item{AIC}{AIC value of the specified model.}

\item{prob}{the estimate of covariates' probability mass function.}

}







\details{
Options \code{d} and \code{K} are necessary.
In the absence of missing data, \code{x} stands for the observed covariates. Refer to Liu et al. (2017) for the EL method.

In the presence of missing data, \code{x} stands for the covariates without missingness, and \code{y} stands for the covariates with missing data (\code{NA}). Refer to Liu et al. (2020+) for the EL method.

Caution should be taken that \code{beta.initial} is NULL if all covariates are considered. When some coefficients of always observed covariates are zero, please specify the corresponding components of \code{beta.initial} as zero and others as none-zero values.
}




\references{
Liu, Y., P. Li, and J. Qin (2017).
Maximum empirical likelihood estimation for abundance in a closed population from capture-recapture data. \emph{Biometrika} \strong{104}(3), 527--543.

Liu, Y., Y. Liu, P. Li, and L. Zhu (2020+). Maximum likelihood abundance estimation from capture-recapture data when covariates are missing at random. Submitted.
}




\examples{
  ## Analysis results of the prinia data set.

  attach(prinia)
  tail.length[tail.length==0] <- NA

  ### EL in the absence of missingness (CC)
  r <- !is.na(tail.length)
  ## Model 1
  (cc1 <- abun.opt(d = number.of.capture, x = cbind(1,fat.index), K = 17) )
  ## Model 2
  (cc2 <- abun.opt(d = number.of.capture, x = cbind(1,wing.index), K = 17) )
  ## Model 3
  (cc3 <- abun.opt(d = number.of.capture, x = cbind(1,fat.index, wing.index), K = 17) )
  ## Model 4
  (cc4 <- abun.opt(d = number.of.capture[r], x = cbind(1, tail.length)[r,], K = 17) )
  ## Model 5
  (cc5 <- abun.opt(d = number.of.capture[r], x = cbind(1,fat.index, tail.length)[r,], K = 17) )
  ## Model 6
  (cc6 <- abun.opt(d = number.of.capture[r], x = cbind(1,wing.index, tail.length)[r,], K = 17) )
  ## Model 7
  (cc7 <- abun.opt(d = number.of.capture[r], x = cbind(1,fat.index, wing.index, tail.length)[r,], K = 17) )

  ### EL in the presence of missingness (EL)
  ## Model 7
  (el7 <- abun.opt(d = number.of.capture, x = cbind(1, fat.index, wing.index), y = tail.length, K=17,
                   beta.initial = c(-10,  1,  0.5, 0.1) ) )
  ## Model 4
  (el4 <- abun.opt(d = number.of.capture, x = cbind(1, fat.index, wing.index), y = tail.length, K = 17,
                   beta.initial = el7$beta*c(1, 0, 0, 1)) )
  ## Model 5
  (el5 <- abun.opt(d = number.of.capture, x = cbind(1, fat.index, wing.index), y = tail.length, K = 17,
                   beta.initial = el7$beta*c(1, 1, 0, 1)) )
  ## Model 6
  (el6 <- abun.opt(d = number.of.capture, x = cbind(1, fat.index, wing.index), y = tail.length, K = 17,
                   beta.initial = el7$beta*c(1, 0, 1, 1)) )

  detach(prinia)

}









