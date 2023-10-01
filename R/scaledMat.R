#' Scale a square matrix
#'
#' Scale a square matrix to have unit diagonal elements.
#' 
#' @param x A square matrix with positive diagonal elements
#'
#' @returns Scaled matrix of x.
#' 
#' @author Min Jin Ha
#' @examples
#'  ###############################
#'  # Simulate data
#'  ###############################
#'  simulation <- simulateData(G = 100, etaA = 0.02, n = 50, r = 10)
#'  dat <- simulation$data[[1L]]
#'  correlation <- scaledMat(x = stats::cov(dat))
#'
#' @importFrom stats cov
#' @export
scaledMat <- function(x) {

  diagx <- diag(x)

  x / sqrt(diagx %o% diagx)
}
