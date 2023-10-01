#' @noRd
#' @keywords internal
svdFunc <- function(x, y, lambda) {

  k <- length(lambda)

  Ds <- tryCatch(svd(x = x), 
                 error = function(e) {
                   stop("unable to obtain SVD\n\t", e$message, call. = FALSE)
                 })

  ld <- length(Ds$d)

  div <- Ds$d*Ds$d + rep(lambda, rep(ld, k))

  rhs <- t(Ds$u) %*% y
  tmp <- drop(Ds$d * rhs) / div
  dim(tmp) <- c(ld, k)

  Ds$v %*% tmp

}
