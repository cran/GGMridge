#' Estimation of empirical null distribution.
#' 
#' Estimation of empirical null distribution using Efron's central matching.
#'
#' @param z A numeric vector of z values following the theoretical normal null 
#'   distribution.
#' @param bins The number of intervals for density estimation of the marginal 
#'   density of z.
#' @param maxQ The maximum degree of the polynomial to be considered for density
#'   estimation of the marginal density of z.
#' @param pct Low and top (pct*100)% tails of z values are excluded to estimate
#'   f(z).
#' @param pct0 Low and top (pct0*100)% tails of z values are excluded to 
#'   estimate f0(z).
#' @param cc The central parts 
#'   \deqn{(\mu - \sigma cc, \mu + \sigma cc)}{(mu - sigma*cc, mu + sigma*cc)}
#'    of the empirical distribution z are used for an estimate of the null 
#'    proportion (eta). 
#' @param plotIt TRUE if density plot is to be produced.
#' 
#' @returns A list containing
#'     \item{correctz }{The corrected z values to follow empirically 
#'                    standard normal distribution.}
#'    \item{correctp }{The corrected p values using the correct z values.}
#'    \item{q }{The chosen degree of polynomial for the estimated 
#'                    marginal density.}
#'    \item{mu0hat }{The location parameter for the normal null 
#'                    distribution.}
#'    \item{sigma0hat }{The scale parameter for the normal null 
#'                    distribution.}
#'    \item{eta }{The estimated null proportion.}
#' @references
#'   Efron, B. 
#'   (2004). 
#'   Large-scale simultaneous hypothesis testing. 
#'   Journal of the American Statistical Association, 99, 96--104.
#'
#'   Ha, M. J. and Sun, W. 
#'   (2014).
#'   Partial correlation matrix estimation using ridge penalty 
#'   followed by thresholding and re-estimation.
#'   Biometrics, 70, 762--770.
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
#'  data <- simulation$data[[1]]
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
#'  w.upper <- which(upper.tri(diag(p)))
#'  partial <- solve(lambda * diag(p) + cor(data))
#'  partial <- (-scaledMat(x = partial))[w.upper]
#'    
#'  ###############################
#'  # get p-values from empirical 
#'  # null distribution 
#'  ###############################
#'  efron.fit <- getEfronp(z = transFisher(x = partial))
#'  
#' @importFrom graphics hist legend lines
#' @importFrom stats coef dnorm glm ks.test lm pnorm quantile
#' @include ksStat.R
#' @export
getEfronp <- function(z, 
                      bins = 120L,
                      maxQ = 9.0,
                      pct = 0.0,
                      pct0 = 0.25,
                      cc = 1.2,
                      plotIt = FALSE) {
  
  stopifnot(
    "`bins` must be > 1" = bins > 1,
    "`pct0` must be >= pct" = pct0 >= pct,
    "`pct must be >= 0" = pct >= 0.0
  )

  pvalues <- 2.0 * {1.0 - stats::pnorm(abs(z))}

  # Z values for marginal density estimation
  v <- stats::quantile(z, probs = c(pct, 1.0 - pct))
  lo <- v[1L]
  hi <- v[2L]

  zz <- pmax(pmin(z, hi), lo)
  
  # X values for marginal density estimation
  breaks <- seq(from = lo, to = hi, length = bins + 1L)
  hcomb <- graphics::hist(x = zz, breaks = breaks, plot = FALSE) 
  counts <- hcomb$counts

  X <- sapply(1L:maxQ, FUN = function(j, x) { x^j }, x = hcomb$mids)
    
  kstest <- rep(NA, maxQ)
  mu0hat <- rep(NA, maxQ)
  sigma0hat <- rep(NA, maxQ)
  eta <- rep(NA, maxQ)

  kstest[1L] <- stats::ks.test(x = pvalues, y = "punif")$statistic
  mu0hat[1L] <- 0.0
  sigma0hat[1L] <- 1.0

  cof <- 0.0
  for (p in 2L:maxQ) {

    fit <- try(stats::glm(counts ~ X[,1L:p], family = "poisson"), 
               silent = TRUE)

    if (inherits(fit, "try-error")) next

    cof <- stats::coef(fit)

    if (any(is.na(cof))) next

    temp <- sapply(0L:p, 
                   FUN = function(i, zz) { zz^i }, 
                   zz = zz)

    fz <- which.max(exp(temp %*% cof))

    mu0hat[p] <- z[fz]

    temp <- sapply(0L:{p-2L}, 
                   FUN = function(i, mu){ {i + 2L} * {i + 1L} * mu^i},
                   mu = mu0hat[p])

    logfzdiff2 <- temp %*% cof[-{1L:2L}]

    if (logfzdiff2 < 0.0) {
      sigma0hat[p] <- {-logfzdiff2}^(-0.5)
      correctZ <- {z - mu0hat[p]} / sigma0hat[p]
      correctp <- 2.0 * {1.0 - stats::pnorm(abs(correctZ))}
      kstest[p] <- ksStat(p = correctp)
    }  
  }

  qq <- which.min(kstest[-1L]) + 1L
  correctZ <- {z - mu0hat[qq]} / sigma0hat[qq]
  correctp <- 2.0 * {1.0 - stats::pnorm(q = abs(correctZ))}
  mu <- mu0hat[qq]
  sigma <- sigma0hat[qq]
  
  # Eta calculation (null proportion)
  fit <- tryCatch(stats::glm(counts ~ X[,1:qq], family = 'poisson'),
                  error = function(e) {
                    stop("glm fit not successful\n\t",
                         e$message, call. = FALSE)
                  })

  logf <- log(fit$fitted.values)

  xmax <- hcomb$mids[which.max(logf)]
  
  v0 <- stats::quantile(z, probs = c(pct0, 1.0 - pct0))
  lo0 <- v0[1L]
  hi0 <- v0[2L]
  idx0 <- which({hcomb$mids > lo0} & {hcomb$mids < hi0})

  if (length(idx0) > 10L) {
    x0 <- hcomb$mids[idx0]
    y0 <- logf[idx0]
    X0 <- cbind(x0 - mu, (x0 - mu)^2)
    fit.null <- tryCatch(stats::lm(y0 ~ X0),
                         error = function(e) {
                           stop("lm fit of null model not successful\n\t",
                                e$message, call. = FALSE)
                         })
    coef.null <- fit.null$coef
  } else {
    idz0 <- which({z > lo0} & {z < hi0})
    z0 <- z[idz0]
    temp <- sapply(0L:qq, FUN = function(i, z) { z^i }, z = z)
    lz <- temp %*% coef(fit)
    y0 <- lz[idz0]
    Z0 <- cbind(z0-mu,(z0-mu)^2)
    fit.null <- tryCatch(stats::lm(y0 ~ Z0), 
                         error = function(e) {
                           stop("lm fit of null model not successful\n\t",
                                e$message, call. = FALSE)
                         })
    coef.null <- fit.null$coef
  }

  # Correct mu and sigma
  mu <- -coef.null[2L] / (2.0 * coef.null[3L]) + mu
  sigma <- sqrt(-0.5 / coef.null[3L])
  lo00 <- mu - cc * sigma
  hi00 <- mu + cc * sigma
  pi <- sum( {z > lo00} & {z < hi00} ) / length(z)
  G0 <- 1.0 - 2.0 * stats::pnorm(lo00, mean = mu, sd = sigma)
  eta <- max(min(pi / G0, 1.0), 0.0)
  
  if (plotIt) {
    oz <- order(zz)
    h <- graphics::hist(x = zz, 
                        breaks = breaks, 
                        xlab = expression(psi), 
                        main = "")
    graphics::lines(x = hcomb$mids, y = fit$fitted.values, col = 3L, lwd = 2.0)
    graphics::lines(x = zz[oz], 
                    y = stats::dnorm(x = zz, mean = mu, sd = sigma)[oz] * 
                      length(zz) * mean(diff(h$breaks)),
                    col = 2L,
                    lty = 2L,
                    lwd = 1.5)
    graphics::legend(x = "topright",
                     legend = c("f", expression(f[0])),
                     col = c(3L,2L),
                     lty = 1L:2L,
                     lwd = c(2.0, 1.5),
                     bty = "n",
                     seg.len = 3L)
  }
  
  list("correctz" = correctZ,
       "correctp" = correctp,
       "q" = qq,
       "mu0hat" = mu,
       "sigma0hat" = sigma,
       "etahat" = eta)
}
