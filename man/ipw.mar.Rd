
\name{ipw.mar}
\alias{ipw.mar}


\title{Inverse-probability-weighting-based (IPW) abundance estimation under the Huggins-Alho (H-A) model}


\description{
This function is used to calculate the IPW estimate and Wald-type confidence interval of the abundance under the H-A model, as well as the estimates of parameters in the H-A model.
}


\usage{

ipw.mar(d, K, x = NULL, y, level =0.95,
  beta.initial = NULL)

}


\arguments{

\item{d}{a vector, the number of times of being captured.}

\item{K}{a number, the number of capture occasions.}

\item{x}{a vector, matrix or data.frame, stands for individual covariates without missingness.}

\item{y}{a vector, matrix or data.frame, stands for individual covariates with missing values (\code{NA}).}

\item{level}{a number, the nominal confidence level (default: 0.95).}

\item{beta.initial}{a vector, with the same length as the parameters in H-A model. 0 means that the coefficient is fixed to zero, and non-zero values are seen as initial values when optimizing.}

}

\value{

A list containing the following components:

\item{n.big}{the IPW estimate of the abundance.}

\item{n.big.se}{Std. Error of \code{n.big}.}

\item{n.big.ci}{the confidence interval of the abundance.}

\item{beta}{the IPW estimates of parameters in the H-A model.}

\item{beta.se}{Std. Error of \code{beta}.}

}







\details{
Options \code{d}, \code{K}, and \code{y} are necessary.

Caution should be taken that \code{beta.initial} is \link{NULL} if all covariates are considered. When some  coefficients of always observed covariates are zero, please specify the corresponding components of \code{beta.initial} as zero and others as none-zero values.
}




\references{
Lee, S.-M., W.-H. Hwang, and J. de Dieu Tapsoba (2016).
Estimation in closed capture-recapture models when covariates are missing at random.
\emph{Biometrics} \strong{72}(4), 1294--1304.
}




\examples{
## Analysis results of the prinia data set.

attach(prinia)

## Model 4
(ipw4 <- ipw.mar( d = number.of.capture, K=17, x = cbind(1, fat.index, wing.index), y = tail.length,
                  beta.initial = c(-10,  0,  0, 0.1) ) )
## Model 5
(ipw5 <- ipw.mar( d = number.of.capture, K=17, x = cbind(1, fat.index, wing.index), y = tail.length,
                  beta.initial = c(-10,  1,  0, 0.1) ) )
## Model 6
(ipw6 <- ipw.mar( d = number.of.capture, K=17, x = cbind(1, fat.index, wing.index), y = tail.length,
                  beta.initial = c(-10,  0,  1, 0.1) ) )
## Model 7
(ipw7 <- ipw.mar( d = number.of.capture, K=17, x = cbind(1, fat.index, wing.index), y = tail.length,
                  beta.initial = c(-10,  1,  0.5, 0.1) ) )

detach(prinia)

}









