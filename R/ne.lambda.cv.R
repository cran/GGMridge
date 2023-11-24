#' Choose the Tuning Parameter of a Ridge Regression Using Cross-Validation
#' 
#' Choose the tuning parameter of a ridge regression using cross-validation.
#' 
#' @param y Length n response vector.
#' @param x n x p matrix for covariates with p variables and n sample size.
#' @param lambda A numeric vector for candidate tuning parameters for a ridge 
#'   regression.
#' @param fold  fold-cross validation used to choose the tuning parameter. 
#'
#' @returns A list containing
#'    \item{lambda }{The selected tuning parameter, which minimizes 
#'      the prediction error. }
#'    \item{spe    }{The prediction error for all of the candidate 
#'      lambda values.}
#'  
#' @references
#'  Ha, M. J. and Sun, W. 
#'  (2014).
#'  Partial correlation matrix estimation using ridge penalty followed 
#'  by thresholding and re-estimation.
#'  Biometrics, 70, 762--770.
#'  
#' @author Min Jin Ha
#' 
#' @examples
#'  p <- 100 # number of variables
#'  n <- 50 # sample size
#'  
#'  ###############################
#'  # Simulate data
#'  ###############################
#'  simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
#'  data <- simulation$data[[1L]]
#'  stddat <- scale(x = data, center = TRUE, scale = TRUE)
#'  
#'  X <- stddat[,-1L,drop = FALSE]
#'  y <- stddat[,1L]
#'  
#'  fit.lambda <- ne.lambda.cv(y = y,
#'                             x = X,
#'                             lambda = seq(from = 0.01, to = 1,by = 0.1),
#'                             fold = 10L)  
#'  
#'  lambda <- fit.lambda$lambda[which.min(fit.lambda$spe)] 
#'
#' @include splitSets.R svdFunc.R
#' @export
ne.lambda.cv <- function(y, x, lambda, fold) {
    
  x <- cbind(y, x)

  n <- nrow(x)

  cv <- {1L:n} %% fold + 1L
  cv <- cbind(cv, sample(1:n))

  k <- length(lambda)
  spe <- numeric(k)

  for (i in 1L:fold) {

    sets <- splitSets(cv = cv, i = i, x = x)

    coef <- svdFunc(x = sets$train[, -1L, drop = FALSE],
                    y = sets$train[, 1L],
                    lambda = lambda)

    spe <- spe + colSums( {sets$test[, 1L] - 
                           sets$test[, -1L, drop = FALSE] %*% coef}^2) 
  }

  list("lambda" = lambda, "spe" = spe / fold)
}
