#' Choose the Tuning Parameter of the Ridge Inverse
#' 
#' Choose the tuning parameter of the ridge inverse by minimizing cross 
#'   validation estimates of the total prediction errors of the p separate 
#'   ridge regressions.
#'
#' @param x An n by p data matrix.
#' @param lambda A numeric vector of candidate tuning parameters.
#' @param fold fold-cross validation is performed.
#'
#' @returns A list containing
#'  \item{lambda }{The selected tuning parameter, which minimizes the 
#'    total prediction errors. }
#'  \item{spe }{The total prediction error for all the candidate 
#'    lambda values.}
#' 
#' @references
#'  Ha, M. J. and Sun, W. 
#'  (2014).
#'  Partial correlation matrix estimation using ridge penalty followed 
#'  by thresholding and re-estimation.
#'  Biometrics, 70, 762--770.
#'
#' @author Min Jin Ha
#' @examples
#'  p <- 100 # number of variables
#'  n <- 50 # sample size
#'  
#'  ###############################
#'  # Simulate data
#'  ###############################
#'  simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
#'  data <- simulation$data[[1L]]
#'  stddata <- scale(x = data, center = TRUE, scale = TRUE)
#'  
#'  ###############################
#'  # estimate ridge parameter
#'  ###############################
#'  lambda.array <- seq(from = 0.1, to = 20, by = 0.1) * (n - 1.0)
#'  fit <- lambda.cv(x = stddata, lambda = lambda.array, fold = 10L)
#'  lambda <- fit$lambda[which.min(fit$spe)] / (n - 1.0)
#'  
#'  ###############################
#'  # calculate partial correlation 
#'  # using ridge inverse
#'  ###############################
#'  partial <- solve(lambda*diag(p) + cor(data))
#'  partial <- -scaledMat(x = partial)
#'
#' @include splitSets.R svdFunc.R
#' @export
lambda.cv <- function(x, lambda, fold) {

  n <- nrow(x)
  p <- ncol(x)

  cv <- {1L:n} %% fold + 1L
  cv <- cbind(cv, sample(1L:n))

  k <- length(lambda)

  spe <- numeric(k)

  for (i in 1L:fold) {

    sets <- splitSets(cv = cv, i = i, x = x)

    speMat <- matrix(0.0, nrow = p, ncol = k)

    for (j in 1L:p) {

      coef <- svdFunc(x = sets$train[,-j,drop = FALSE],
                      y = sets$train[,j],
                      lambda = lambda)

      speMat[j,] <- colSums( {sets$test[, j] - 
                              sets$test[, -j, drop = FALSE] %*% coef}^2 )
    }

    spe <- spe + colSums(speMat)

  }

  list("lambda" = lambda, "spe" = spe / {fold * p})
}
