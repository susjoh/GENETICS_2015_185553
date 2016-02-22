pin <- function (object, transform) 
{
  pframe <- as.list(object$gammas)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
                 pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  tname <- if (length(transform) == 3) {
    transform[[2]]
  } else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}
