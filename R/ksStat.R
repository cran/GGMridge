#'The Kolmogorov-Smirnov Statistic for p-Values
#'
#' Calculates the Kolmogorov-Smirnov statistic for p-values
#' 
#' @param p A numeric vector with p-values.
#'
#' @returns Kolmogorov-Smirnov statistic
#' @author Min Jin Ha
#' @importFrom stats punif runif
#' @examples
#'  p <- stats::runif(100)
#'  ksStat(p = p)
#'  ks.test(p, y = "punif") # compare with ks.test
#'  
#' @export
ksStat <- function(p) {

  p <- p[!is.na(p)]
  n <- length(p)

  p <- stats::punif(q = sort(p)) - {0L:{n - 1L}} / n

  max(c(p, 1.0/n - p))

}
