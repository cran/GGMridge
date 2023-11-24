#' Estimation of the mixture distribution using EM algorithm
#' 
#' Estimation of the parameters, null proportion, and 
#'  degrees of freedom of the exact null density in the 
#'  mixture distribution.
#'
#'
#' @param p A numeric vector representing partial correlation coefficients.
#' @param eta0 An initial value for the null proportion; 1-eta0 is the non-null 
#'   proportion.
#' @param df  An initial value for the degrees of freedom of the exact null 
#'   density.
#' @param tol The tolerance level for convergence.
#'
#' @returns A list object containing
#'     \item{df   }{Estimated degrees of freedom of the null density.}
#'     \item{eta0 }{Estimated null proportion.}
#'     \item{iter }{The number of iterations required to reach convergence.}
#' @references
#'  Schafer, J. and Strimmer, K. 
#'  (2005). 
#'  An empirical Bayes approach to inferring large-scale gene association 
#'  networks. 
#'  Bioinformatics, 21, 754--764.
#'  
#' @author Min Jin Ha
#'
#' @importFrom stats dunif
#' @export
EM.mixture <- function(p, eta0, df, tol) {

  f0 <- function(x, df) {
          (1.0 - x^2)^{0.5*{df - 3.0}} * 
          (1.0 / beta(a = 0.5, b = (df - 1.0)*0.5))
        }
  fa <- function(x){ stats::dunif(x = x, min = -1.0, max = 1.0) }
  
  E <- length(p)
  eps <- max(1000.0, tol*1.1)
  i <- 0L
  while( eps > tol ) {
    i <- i + 1L
    temp <- eta0 * f0(x = p, df = df)

    # conditional expection of the missing
    Ak <- temp / {temp + {1.0 - eta0} * fa(x = p)} 

    sumAk <- sum(Ak)
    
    # condition expection of first derivative of loglikelihood
    expect1 <- {sum(log(1.0 - p^2) * Ak) + 
                {digamma(x = 0.5*df) - digamma(x = 0.5*(df-1.0))} * 
                 sumAk} * 0.5
    
    # condition expection of second derivative of loglikelihood
    expect2 <- {trigamma(x = 0.5*df) - trigamma(x = 0.5*(df-1.0))} * 
               sumAk  * 0.25
    tempEta0 <- sumAk / E 

    # using newton-raphson
    tempDF <- df - {1.0 / expect2} * expect1
    eps <- max( abs(eta0 - tempEta0), abs(df - tempDF) )
    eta0 <- tempEta0
    df <- tempDF
  }

  list("df" = df, 
       "eta0" = eta0, 
       "iter" = i)
}
