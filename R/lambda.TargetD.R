#' Shrinkage Estimation of a Covariance Matrix Toward an Identity Matrix 
#' 
#' Estimation of a weighted average of a sample covariance (correlation) matrix 
#' and an identity matrix.
#' 
#' An analytical approach to the estimate ridge parameter.
#' 
#' @param x Centered data for covariance shrinkage and standardized data for 
#'   correlation shrinkage.
#'
#' @returns The estimates of shrinkage intensity.
#'
#' @references
#' Schafer, J. and Strimmer, K.
#'  (2005). 
#'  A shrinkage approach to large-scale covariance matrix estimation 
#'  and implications for functional genomics. 
#'  Statistical Applications in Genetics and Molecular Biology, 4, 32.
#'  
#'  Ha, M. J. and Sun, W. 
#'  (2014).
#'  Partial correlation matrix estimation using ridge penalty followed 
#'  by thresholding and re-estimation.
#'  Biometrics, 70, 762--770.
#'
#' @author Min Jin Ha
#' 
#' @examples
#' ###############################
#'  # Simulate data
#'  ###############################
#'  simulation <- simulateData(G = 100, etaA = 0.02, n = 50, r = 10)
#'  dat <- simulation$data[[1L]]
#'  stddat <- scale(x = dat, center = TRUE, scale = TRUE)
#'  
#'  shrinkage.lambda <- lambda.TargetD(x = stddat)
#'  
#'  ###############################
#'  # the ridge parameter
#'  ###############################
#'  ridge.lambda <- shrinkage.lambda / (1.0 - shrinkage.lambda)
#'  
#'  ###############################
#'  # partial correlation matrix
#'  ###############################
#'  partial <- solve(cor(dat) + ridge.lambda * diag(ncol(dat)))
#'  partial <- -scaledMat(x = partial)
#'  
#' @export
lambda.TargetD <- function(x) {

  n <- nrow(x)

  wmean2 <- {t(x) %*% x}^2

  x2 <- x^2
  x2 <- t(x2) %*% x2

  diag.w <- which(row(wmean2) == col(wmean2))

  varw <- sum(({ x2  - wmean2 / n} * n / {{n-1}^3})[-diag.w])

  S2 <- sum(wmean2[-diag.w] / {n - 1L}^2)

  max(0.0, min( varw / S2, 1.0))
}
