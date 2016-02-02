\name{EM.mixture}
\alias{EM.mixture}
\title{ Estimation of the mixture distribution using EM algorithm }
\description{ Estimation of the parameters, null proportion, and 
              degrees of freedom of the exact null density in the 
              mixture distribution.}
\usage{ EM.mixture(p, eta0, df, tol) }
\arguments{
  \item{p    }{A numeric vector representing partial correlation 
               coefficients.}
  \item{eta0 }{An initial value for the null proportion; 1-eta0 is 
               the non-null proportion.}
  \item{df   }{An initial value for the degrees of freedom of the 
               exact null density.}
  \item{tol  }{The tolerance level for convergence.}
}
\value{
  \item{df   }{Estimated degrees of freedom of the null density.}
  \item{eta0 }{Estimated null proportion.}
  \item{iter }{The number of iterations required to reach convergence.}
}
\references{
Schafer, J. and Strimmer, K. 
(2005). 
An empirical Bayes approach to inferring large-scale gene 
association networks. 
Bioinformatics, 21, 754--764.
}
\author{ Min Jin Ha }
