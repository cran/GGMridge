\name{ksStat}
\alias{ksStat}
\title{ The Kolmogorov-Smirnov Statistic for p-Values }
\description{ Calculates the Kolmogorov-Smirnov statistic for p-values}
\usage{ ksStat(p) }
\arguments{
  \item{p}{A numeric vector with p-values.}
}
\value{ Kolmogorov-Smirnov statistic }
\author{ Min Jin Ha }
\examples{
  p <- runif(100)
  ksStat(p = p)
  ks.test(p, y = "punif") # compare with ks.test
}
