\name{getEfronp}
\alias{getEfronp}
\title{ Estimation of empirical null distribution. }
\description{ Estimation of empirical null distribution using Efron's 
              central matching. }
\usage{ getEfronp(z, bins = 120L, maxQ = 9, pct = 0.0, pct0 = 0.25, 
        cc = 1.2, plotIt = FALSE) }
\arguments{
  \item{z     }{A numeric vector of z values following the theoretical 
                normal null distribution.}
  \item{bins  }{The number of intervals for density estimation of the 
                marginal density of z.}
  \item{maxQ  }{The maximum degree of the polynomial to be considered 
                for density estimation of the marginal density of z.}
  \item{pct   }{Low and top (pct*100)% tails of z values are excluded 
                to estimate f(z).}
  \item{pct0  }{Low and top (pct0*100)% tails of z values are 
                excluded to estimate f0(z).}
   \item{cc   }{The central parts (mu - sigma*cc, mu + sigma*cc) of 
                the empirical distribution z are used for an 
                estimate of the null proportion (eta). }
  \item{plotIt}{TRUE if density plot is to be produced.}
}
\value{
  \item{correctz  }{The corrected z values to follow empirically 
                    standard normal distribution.}
  \item{correctp  }{The corrected p values using the correct z values.}
  \item{q         }{The chosen degree of polynomial for the estimated 
                    marginal density.}
  \item{mu0hat    }{The location parameter for the normal null 
                    distribution.}
  \item{sigma0hat }{The scale parameter for the normal null 
                    distribution.}
  \item{eta       }{The estimated null proportion.}

}
\references{
Efron, B. 
(2004). 
Large-scale simultaneous hypothesis testing. 
Journal of the American Statistical Association, 99, 96--104.

Ha, M. J. and Sun, W. 
(2014).
Partial correlation matrix estimation using ridge penalty 
followed by thresholding and re-estimation.
Biometrics, 70, 762--770.
}
\author{ Min Jin Ha }
\examples{ 
  p <- 100 # number of variables
  n <- 50 # sample size
	
  ###############################
  # Simulate data
  ###############################
  simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
  data <- simulation$data[[1]]
  stddata <- scale(x = data, center = TRUE, scale = TRUE)
    
  ###############################
  # estimate ridge parameter
  ###############################
  lambda.array <- seq(from = 0.1, to = 20, by = 0.1) * (n - 1.0)
  fit <- lambda.cv(x = stddata, lambda = lambda.array, fold = 10L)
  lambda <- fit$lambda[which.min(fit$spe)] / (n - 1.0)
    
  ###############################
  # calculate partial correlation
  # using ridge inverse
  ###############################
  w.upper <- which(upper.tri(diag(p)))
  partial <- solve(lambda * diag(p) + cor(data))
  partial <- (-scaledMat(x = partial))[w.upper]
    
  ###############################
  # get p-values from empirical 
  # null distribution 
  ###############################
  efron.fit <- getEfronp(z = transFisher(x = partial))
}

