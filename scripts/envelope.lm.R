## ID: envelope.lm.R, extracted from Venables and Ripley (2000)
## S Programming, Section 7.3

envelope <- function(object, nsamples = 1000, alpha = 0.05) {
   n <- length(r <- resid(object))
   p <- length(coef(object))
   Y <- scale(qr.resid(qr(model.matrix(object)),
   		matrix(rnorm(n * nsamples), n, nsamples)),
   		F, T) * sqrt((n - 1) / (n - p))
   Y[,] <- Y[order(col(Y),Y)]
   Y <- matrix(Y[order(row(Y),Y)], nsamples, n)
   x0 <- quantile(1:nsamples, c(alpha / 2, 1 - alpha / 2))
   if (all(x0 %% 1 == 0)) elim <- t(Y[x0,])
   else {
      x1 <- c(floor(x0), ceiling(x0))
      elim <- cbind(Y[x1[1],] + (Y[x1[3],] - Y[x1[1],]) /
      		(x1[3] - x1[1]) * (x0[1] - x1[1]),
      		Y[x1[2],] + (Y[x1[4],] - Y[x1[2],]) /
      		(x1[4] - x1[2]) * (x0[2] - x1[2]))
   }
   res <- sort(r)
   res <- res / sqrt(sum(res^2) / (n - p))
   res <- structure(
   		list(res = res, elim = elim),
   		label = deparse(object$call$formula))
   class(res) <- "envelope"
   res
}

plot.envelope <- function(x, ...) {
   n <- length(x$res)
   ylim <- range(x$res[c(1,n)], x$elim[c(1,n),])
   nscores <- qnorm(ppoints(n))
   oldpar <- par(pty = "s"); on.exit(par(oldpar))
   plot(nscores, x$res, pch = 1, ylim = ylim,
   	xlab = "Normal scores", ylab = "Sorted residuals",
   	main = attr(x, "label"))
   lines(nscores, x$elim[,1]); lines(nscores, x$elim[,2])
   invisible(x)
}

print.envelope <- function(x, ...) {
   lo <- x$res < x$elim[,1]
   hi <- x$res > x$elim[,2]
   flash <- rep("", length(lo))
   flash[lo] <- "<"; flash[hi] <- ">"
   if (any(lo | hi))
      print(cbind(do.call("data.frame", x), flash = flash)[lo | hi,])
   else cat("All points within envelope\n")
   invisible(x)
}
