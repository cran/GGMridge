% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne.lambda.cv.R
\name{ne.lambda.cv}
\alias{ne.lambda.cv}
\title{Choose the Tuning Parameter of a Ridge Regression Using Cross-Validation}
\usage{
ne.lambda.cv(y, x, lambda, fold)
}
\arguments{
\item{y}{Length n response vector.}

\item{x}{n x p matrix for covariates with p variables and n sample size.}

\item{lambda}{A numeric vector for candidate tuning parameters for a ridge 
regression.}

\item{fold}{fold-cross validation used to choose the tuning parameter.}
}
\value{
A list containing
   \item{lambda }{The selected tuning parameter, which minimizes 
     the prediction error. }
   \item{spe    }{The prediction error for all of the candidate 
     lambda values.}
}
\description{
Choose the tuning parameter of a ridge regression using cross-validation.
}
\examples{
 p <- 100 # number of variables
 n <- 50 # sample size
 
 ###############################
 # Simulate data
 ###############################
 simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
 data <- simulation$data[[1L]]
 stddat <- scale(x = data, center = TRUE, scale = TRUE)
 
 X <- stddat[,-1L,drop = FALSE]
 y <- stddat[,1L]
 
 fit.lambda <- ne.lambda.cv(y = y,
                            x = X,
                            lambda = seq(from = 0.01, to = 1,by = 0.1),
                            fold = 10L)  
 
 lambda <- fit.lambda$lambda[which.min(fit.lambda$spe)] 

}
\references{
Ha, M. J. and Sun, W. 
 (2014).
 Partial correlation matrix estimation using ridge penalty followed 
 by thresholding and re-estimation.
 Biometrics, 70, 762--770.
}
\author{
Min Jin Ha
}
