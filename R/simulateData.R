#' Generate Simulation Data from a Random Network.
#' 
#' Generate a random network where both the network structure and the partial 
#'   correlation coefficients are random. The data matrices are generated from 
#'   multivariate normal distribution with the covariance matrix corresponding 
#'   to the network.
#'   
#' @param G The number of variables (vertices).
#' @param etaA The proportion of non-null edges among all the G(G-1)/2 edges.
#' @param n The sample size.
#' @param r The number of replicated G by N data matrices.
#' @param dist A function which indicates the distribution of sample. 
#'    "mvnorm" is multivariate normal distribution and 
#'    "mvt" is multivariate t distribution with df=2. 
#'    The default is set by "mvnorm".
#'    
#' @returns A list containing
#'   \itemize{
#'     \item{data }{ a list, each element containing an n X G matrix of 
#'     simulated data.}
#'     \item{true.partialcor }{ The partial correlation matrix which the 
#'     datasets are generated from.}
#'     \item{truecor.scaled }{ The covariance matrix calculted from the 
#'     partial correlation matrix.}
#'     \item{sig.node }{ The indices of nonzero upper triangle 
#'     elements of partial correlation matrix.}
#'  }
#'  
#' @references
#'  Schafer, J. and Strimmer, K. 
#'  (2005). 
#'  An empirical Bayes approach to inferring large-scale gene 
#'  association networks. 
#'  Bioinformatics, 21, 754--764.
#'  
#' @author Min Jin Ha
#' @examples
#'  simulation <- simulateData(G = 100, etaA = 0.02, n = 50, r = 10)
#'
#' @importFrom mvtnorm rmvnorm rmvt
#' @importFrom stats runif
#' @export
simulateData <- function(G, etaA, n, r, dist = "mvnorm") {

  expression <- list()
  dist <- tolower(dist)

  partialcov <- matrix(0.0, nrow = G, ncol = G)

  tri.w <- which(upper.tri(partialcov))

  no.sig <-  ceiling(etaA * length(tri.w))

  sig.node <- sample(tri.w, no.sig)

  partialcov[sig.node] <- stats::runif(n = no.sig, min = -1.0, max = 1.0)

  partialcov <- t(partialcov)  + partialcov

  diag(partialcov) <- colSums(abs(partialcov)) + 0.0001

  dpc <- diag(partialcov)
  partialcor <- partialcov / sqrt(dpc %o% dpc)
  invcor <- - partialcor
  diag(invcor) <- 1.0
  truecor <- tryCatch(solve(invcor), 
                      error = function(e) {
                        stop("unable to invert correlation matrix\n\t",
                             e$message, call. = FALSE)
                      })

  dtc <- diag(truecor)
  truecor.scaled <- truecor/ sqrt(dtc %o% dtc)

  if( dist == "mvnorm" ) {
    for( i in 1L:r ) {
      expression[[i]] <- mvtnorm::rmvnorm(n = n, mean = rep(0,G), truecor)
    }
  } else if (dist == "mvt") {
    for( i in 1L:r ) {
      expression[[i]] <- mvtnorm::rmvt(n = n, sigma = truecor, df = 2L)
    }
  }

  list("data" = expression,
       "true.partialcor" = partialcor,
       "truecor.scaled" = truecor.scaled,
       "sig.node" = sig.node)
}

