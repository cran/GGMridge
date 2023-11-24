#' Estimation of Partial Correlation Matrix Using p Separate Ridge Regressions.
#' 
#' The partial correlation matrix is estimated by p separate ridge regressions 
#'   with the parameters selected by cross validation.
#'   
#' @param x n x p data matrix; n is the # of samples and p is the # of variables.
#' @param fold Ridge parameters are selected by fold-cross validations 
#'   separately for each regression.
#' @param lambda The candidate ridge parameters for all p ridge regressions.
#' @param verbose TRUE/FALSE; if TRUE, print the procedure.
#'
#' @returns A list containing
#'    \item{R }{ The partial correlation matrix.}
#'    \item{lambda.sel }{ The selected tuning parameters for p ridge regressions.}
#'  
#' @references
#'    Ha, M. J. and Sun, W. 
#'    (2014).
#'    Partial correlation matrix estimation using ridge penalty followed 
#'    by thresholding and re-estimation.
#'    Biometrics, 70, 762--770.
#'    
#' @author Min Jin Ha
#' @examples
#'    p <- 100 # number of variables
#'    n <- 50 # sample size
#'    
#'    ###############################
#'    # Simulate data
#'    ###############################
#'    simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
#'    data <- simulation$data[[1L]]
#'    stddata <- scale(x = data, center = TRUE, scale = FALSE)
#'    
#'    ###############################
#'    # estimate ridge parameter
#'    ###############################
#'    w.upper <- which(upper.tri(diag(p)))
#'    
#'    lambda.array <- seq(from = 0.1, to = 20, by=0.1) * (n-1.0)
#'    partial.sep <-  R.separate.ridge(x = stddata,
#'                                     lambda = lambda.array,
#'                                     fold = 5L,
#'                                     verbose = TRUE)$R[w.upper]
#'  
#' @importFrom MASS lm.ridge
#' @importFrom stats coef
#' @include ne.lambda.cv.R scaledMat.R                                
#' @export
R.separate.ridge <- function(x, fold, lambda, verbose = FALSE) {
  
  stopifnot(
    "`x` must be a matrix" = !missing(x) && is.matrix(x),
    "more than 1 sample is required in `x`" = nrow(x) > 1L,
    "more than 1 covariate is required in `x`" = ncol(x) > 1L
  )
    
  n <- nrow(x)
  p <- ncol(x)
  x <- scale(x = x, center = TRUE, scale = TRUE)
    
  coefNI <- matrix(0.0, nrow = p, ncol = p)
  diag.coefNI <- numeric(p)
  lambda.sel <- numeric(p)

  for (i in 1L:p) {

    if (verbose) message("variable = ", i, " ", date())

    tempX <- x[, -i, drop = FALSE]
    tempY <- x[, i]

    fit.lambda <- ne.lambda.cv(y = tempY,
                               x = tempX,
                               lambda = lambda,
                               fold = fold) 

    lambdai <- fit.lambda$lambda[which.min(fit.lambda$spe)]
    lambda.sel[i] <- lambdai

    ridgefit <- tryCatch(MASS::lm.ridge(tempY ~ tempX - 1, lambda = lambdai), 
                         error = function(e) {
                           stop("ridge regression not successful\n\t",
                                e$message, call. = FALSE)
                         })

    predY <- tempX %*% stats::coef(ridgefit) 
    Ds <- tryCatch(svd(tempX),
                   error = function(e) {
                     stop("unable to obtain SVD\n\t", e$message, call. = FALSE)
                   })
    dvals2 <- Ds$d * Ds$d

    diag.coefNI[i] <- {n - sum(dvals2 / {dvals2 + lambdai})} / 
                      sum({tempY - predY}^2)
    coefNI[i,-i] <- -diag.coefNI[i] * stats::coef(ridgefit)
  }

  absc <- abs(coefNI)
  tmp1 <- sqrt(absc * t(absc))
  tmp2 <- sign(coefNI) * upper.tri(coefNI)
  tmp3 <- (tmp2 + t(tmp2)) * tmp1
  diag(tmp3) <- diag.coefNI
  R <- - scaledMat(x = tmp3)
  diag(R) <- 1.0

  list("R" = R, "lambda.sel" = lambda.sel)
}
