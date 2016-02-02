\name{scaledMat}
\alias{scaledMat}
\title{Scale a square matrix}
\description{ Scale a square matrix to have unit diagonal elements.}
\usage{ scaledMat(x) }
\arguments{
  \item{x}{ A square matrix with positive diagonal elements}
}
\value{ Scaled matrix of x }
\author{ Min Jin Ha }
\examples{
  ###############################
  # Simulate data
  ###############################
  simulation <- simulateData(G = 100, etaA = 0.02, n = 50, r = 10)
  dat <- simulation$data[[1L]]
  correlation <- scaledMat(x = cov(dat))
}
