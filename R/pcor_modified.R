pcor.modified <- function (z, 
                           method = c("pearson", "kendall", "spearman"), 
                           use = c("complete.obs","pairwise.complete.obs")) 
{
  method <- match.arg(method)
  use <- match.arg(use)
  if (is.data.frame(z)) 
    z <- as.matrix(z)
  if (!is.matrix(z)) 
    stop("supply a matrix-like 'z'")
  if (!(is.numeric(z) || is.logical(z))) 
    stop("'z' must be numeric")
  stopifnot(is.atomic(z))
  n <- dim(z)[1]
  gp <- dim(z)[2] - 2
  cvz <- cov(z, method = method, use = use)
  icvz <- solve(cvz)
  pcor <- -cov2cor(icvz)
  diag(pcor) <- 1
  if (method == "kendall") {
    statistic <- pcor/sqrt(2 * (2 * (n - gp) + 5)/(9 * (n - 
                                                          gp) * (n - 1 - gp)))
    p.value <- 2 * pnorm(-abs(statistic))
  }
  else {
    statistic <- pcor * sqrt((n - 2 - gp)/(1 - pcor^2))
    p.value <- 2 * pnorm(-abs(statistic))
  }
  diag(statistic) <- 0
  diag(p.value) <- 0
  list(estimate = pcor, p.value = p.value, statistic = statistic, 
       n = n, gp = gp, method = method)
}