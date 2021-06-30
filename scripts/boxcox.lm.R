## ID: boxcox.lm.R, last updated 2021-06-29, F.Osorio

boxcox.lm <- function(x, y, lambda)
{ # Box-Cox estimation procedure
  boxcox <- function(y, lambda) {
    # transformación Box-Cox
    n <- length(y)
    lambda <- rep(lambda, n) # a little trick
    z <- ifelse(lambda != 0., (y^lambda - 1.) / lambda, log(y))
    z
  }

  n <- nrow(x)
  p <- ncol(x)

  k <- length(lambda)
  RSS <- rep(0, k)
  logLik <- rep(0, k)
  for (i in 1:k) {
    geom <- geomean(y) # nice function from fastmatrix
    z <- boxcox(y, lambda = lambda[i])
    z <- z / geom^(lambda[i] - 1.)
    fm <- ols.fit(x, z, method = "sweep") # 'x' already have an intercept
    RSS[i] <- fm$RSS
    logLik[i] <- -.5 * n * log(2 * pi) - .5 * n * log(RSS[i] / n) - .5 * n
  }
  idx <- order(RSS)[1]
  opt <- lambda[idx]

  obj <- list(lambda = lambda, RSS = RSS, logLik = logLik, opt = opt)
  obj
}
