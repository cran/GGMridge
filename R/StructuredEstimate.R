#' Estimation of Partial Correlation Matrix Given Zero Structure.
#' 
#' Estimation of nonzero entries of the partial correlation matrix given zero 
#'   structure.
#'   
#' @param x n by p data matrix with the number of variables p and sample size n.
#' @param E The row and column indices of zero entries of the partial 
#'   correlation matrix.
#'
#' @returns A list containing
#'     \item{R }{The partial correlation matrix.}
#'     \item{K }{The inverse covariance matrix.}
#'     \item{RSS }{The residual sum of squares.}
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
#'  lambda.array <- seq(from = 0.1, to = 20, by = 0.1) * (n-1.0)
#'  fit <- lambda.cv(x = stddata, lambda = lambda.array, fold = 10L)
#'  lambda <- fit$lambda[which.min(fit$spe)]/(n-1)
#'  
#'  ###############################
#'  # calculate partial correlation 
#'  # using ridge inverse
#'  ###############################
#'  w.upper <- which(upper.tri(diag(p)))
#'  
#'  partial <- solve(lambda * diag(p) + cor(data))
#'  partial <- (-scaledMat(x = partial))[w.upper]
#'  
#'  ###############################
#'  # get p-values from empirical 
#'  # null distribution 
#'  ###############################
#'  efron.fit <- getEfronp(z = transFisher(x = partial), 
#'                         bins = 50L, 
#'                         maxQ = 13)
#'  
#'  ###############################
#'  # estimate the edge set of 
#'  # partial correlation graph with 
#'  # FDR control at level 0.01
#'  ###############################
#'  w.array <- which(upper.tri(diag(p)),arr.ind=TRUE)
#'  th <- 0.01
#'  wsig <- which(p.adjust(efron.fit$correctp, method="BH") < th )
#'  E <- w.array[wsig,]
#'  dim(E)
#'  
#'  ###############################
#'  # structured estimation
#'  ###############################
#'  fit <- structuredEstimate(x = stddata, E = E)
#'  th.partial <- fit$R
#'
#' @importFrom stats var
#' @include scaledMat.R
#' @export
structuredEstimate <- function(x, E) {

  E <- matrix(E, ncol = 2L)

  n <- nrow(x)
  p <- ncol(x)

  RSS <- numeric(p)

  K <- matrix(0.0, nrow = p, ncol = p)
  diag(K) <- apply(x, 2L, FUN = stats::var)

  for (i in 1L:p) {

    y <- x[, i] 

    ne <- c( E[{E[, 1L] == i}, 2L], E[{E[, 2L] == i}, 1L])
    lne <- length(ne)

    if (lne == 0L) {
      RSS[i] <- y %*% y
    } else {
      D <- x[, ne, drop = FALSE] 
      svdD <- svd(x = D)
      coef <- svdD$v %*% {t(svdD$v) / svdD$d^2} %*% t(D) %*% y
      RSS[i] <- sum({y - D %*% coef}^2)
      resid.var <- {n - lne} / RSS[i]
      coef.cl <- c(resid.var, - coef * resid.var)
      id.cl <- c(i, ne)  
      o.cl <- order(id.cl)
      K[i,id.cl[o.cl]] <- coef.cl[o.cl]
    }
  }

  absK <- abs(K)
  tmp1 <- sqrt(absK * t(absK))
  tmp2 <- sign(K) * upper.tri(K)
  tmp3 <- {tmp2 + t(tmp2)} * tmp1
  diag(tmp3) <- diag(K)
  R <- - scaledMat(x = tmp3)
  diag(R) <- 1.0

  list("R" = R, "K" = tmp3, "RSS" = RSS)
}
