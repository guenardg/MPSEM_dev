## AUxiliary functions

OUEvolve <- function(x, from, opt, alpha, sigma) {
  n <- length(x)
  y <- numeric(n)
  y[1L] <- from
  for(i in 2L:n) {
    l <- x[i] - x[i - 1L]
    w <- if(alpha) exp(-alpha*l) else 1
    s <- if(alpha) sigma*sqrt((1 - exp(-2*alpha*l))/(2*alpha)) else sigma*sqrt(l)
    y[i] <- rnorm(1L, w*y[i - 1L] + opt*(1 - w), s)
  }
  y
}
