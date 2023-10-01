#' Choose the Tuning Parameter of the Ridge Inverse and 
#'  Thresholding Level of the Empirical p-Values.
#'  
#' Calculate total prediction error for test data after fitting partial 
#'   correlations from train data for all values of lambda and pcut.
#'
#' @param train An n x p data matrix from which the model is fitted.
#' @param test An m x p data matrix from which the model is evaluated.
#' @param lambda A vector of candidate tuning parameters.
#' @param pcut A vector of candidate cutoffs of pvalues.
#'
#' @returns Total prediction error for all the candidate lambda and pvalue 
#'   cutoff values.
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
#'  
#'  ###############################
#'  # Split into train/test sets
#'  ###############################
#'  testindex <- sample(1L:n, 10L)
#'  
#'  train <- data[-testindex,,drop = FALSE]
#'  stdTrain <- scale(x = train, center = TRUE, scale = TRUE)
#'  
#'  test <- data[testindex,,drop = FALSE]
#'  stdTest <- scale(x = test, center = TRUE, scale = TRUE)
#'  
#'  ###############################
#'  # Calculate total prediction 
#'  # errors for all candidate 
#'  # lambda and p-value cutoffs
#'  ###############################
#'  lambda.array <- seq(from = 0.1, to = 5, length = 10) * (n - 1.0)
#'  pcut.array <- seq(from = 0.01, to = 0.05, by = 0.01)
#'  tpe <- lambda.pcut.cv1(train = stdTrain,
#'                         test = stdTest,
#'                         lambda = lambda.array,
#'                         pcut = pcut.array)
#'   
#' @include getEfronp.R scaledMat.R StructuredEstimate.R transFisher.R
#' @importFrom stats cor                      
#' @export
lambda.pcut.cv1 <- function(train, test, lambda, pcut) {

  p <- ncol(test)
  dp <- diag(p)

  k <- length(lambda)
  lpc <- length(pcut)

  w.upper <- which(upper.tri(dp))
  w.array <- which(upper.tri(dp), arr.ind = TRUE)
  lw <- length(w.upper)

  tunings <- matrix(NA, ncol = 2L, nrow = k * lpc)

  S <- stats::cor(train)

  R <- matrix(NA, nrow = lw, ncol = k)

  kk <- 0L
  for (la in lambda) {

    kk <- kk + 1L

    tmp <- tryCatch(solve(S + la * dp), 
                    error = function(e) {
                      stop("unable to invert matrix\n\t", 
                           e$message, call. = FALSE)
                    })
    
    R[, kk] <- -scaledMat(x = tmp)[w.upper]

  }

  transR <- transFisher(x = R)

  efron <- sapply(1L:k,
                  FUN = function(i) { getEfronp(z = transR[,i])$correctp})

  risk <- matrix(NA, nrow = k, ncol = lpc, 
                 dimnames = list(lambda, pcut))

  for (i in 1L:lpc) {

    thR <- R * {efron < pcut[i]}

    coef.mat <- matrix(NA, nrow = p * p, ncol = k)

    for (kk in 1L:k) {

      notzero <- abs(thR[,kk]) > 1.5e-8

      E <- w.array[notzero]
      fit <- structuredEstimate(x = train, E = E)     
      temp <- sqrt(diag(fit$K))
      coef.mat[,kk] <- diag(temp) %*% fit$R %*% diag(1.0/temp)
    }
   
    CV <- 0.0
    for (j in 1L:p) {
      jp <- j * p
      y <- test[, j]
      D <- test[, -j, drop = FALSE]
      loss <- sapply(1L:k, 
                     FUN = function(i, y, D, coef) {
                             (y - D %*% (coef[,i])[-1L])^2
                           },
                     y = y,
                     D = D,
                     coef = coef.mat[{jp - p + 1L}:jp, , drop = FALSE])

      CV <- CV + colSums(loss) 
    }

    risk[,i] <- CV / length(test) 
  }

  risk
}
