#' Choose the Tuning Parameter of the Ridge Inverse and 
#'   Thresholding Level of the Empirical p-Values
#'
#' Choose the tuning parameter of the ridge inverse and p-value cutoff by 
#'   minimizing cross validation estimates of the total prediction errors of 
#'   the p separate ridge regressions.
#'
#' @param x n by p data matrix.
#' @param lambda A vector of candidate tuning parameters.
#' @param pcut A vector of candidate cutoffs of pvalues.
#' @param fold fold-cross validation is performed.
#'
#' @returns The total prediction errors for all lambda (row-wise) and pcut 
#'   (column-wise)
#'
#' @references
#'  Ha, M. J. and Sun, W. 
#'  (2014).
#'  Partial correlation matrix estimation using ridge penalty followed by 
#'  thresholding and re-estimation.
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
#'  stddata <- scale(x = data, center = TRUE, scale = TRUE)
#'  
#'  ###############################
#'  # Selection of a lambda and a 
#'  # p-value cutoff
#'  ###############################
#'  lambda.array <- seq(from = 0.1, to = 5, length = 10) * (n-1.0)
#'  pcut.array <- seq(from = 0.01, to = 0.05, by = 0.01)
#'  tpe <- lambda.pcut.cv(x = stddata,
#'                        lambda = lambda.array,
#'                        pcut = pcut.array,
#'                        fold = 3L)
#'  w.mintpe <- which(tpe == min(tpe), arr.ind = TRUE)
#'  lambda <- lambda.array[w.mintpe[1L]]
#'  alpha <- pcut.array[w.mintpe[2L]]
#'  
#' @include lambda.pcut.cv1.R splitSets.R
#' @export
lambda.pcut.cv <- function(x, lambda, pcut, fold = 10L) {

  n <- nrow(x)

  cv <- 1L:n %% fold + 1L
  cv <- cbind(cv, sample(1:n)) 

  PE <- matrix(0.0,
               nrow = length(lambda),
               ncol = length(pcut),
               dimnames = list(lambda, pcut))

  k <- length(lambda)

  message("cv: ", appendLF = FALSE)
  for (i in 1L:fold) {
    message(i, " ", appendLF = FALSE)

    sets <- splitSets(cv = cv, i = i, x = x)

    PE <- PE + lambda.pcut.cv1(train = sets$train,
                               test = sets$test,
                               lambda = lambda, 
                               pcut = pcut)
  }
  message("Done!")
  
  PE / fold
}
