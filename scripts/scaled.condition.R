## ID: scaled.condition.R, last updated 2021-06-13, F.Osorio

scaled.condition <- function(x)
{ # scaled condition number
  colScales <- apply(x, 2, function(x) sum(x^2))
  z <- scale(x, center = FALSE, scale = sqrt(colScales))
  d <- svd(z)$d
  p <- length(d)
  cn <- d[1] / d[p]
  obj <- list(condition = cn, values = d, x.scaled = z)
  obj
}
