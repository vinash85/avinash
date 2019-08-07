
#' calc.cor calulate correlation and approximate P-value of two matrix
#'
#' @param x : matrix
#' @param y : matrix
#' @param method : a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param method : matrix
#' @param use : "an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#'
#' @return estimate=correlation coefficient, P=P-value
#' @export
#'
#' @examples
calc.cor = function(x, y, method = "spearman", use = "pairwise.complete.obs"){
  h=suppressWarnings(cor(x, y, method = method, use = use))
  aa = (!is.na(x)) + 0
  bb = (!is.na(y)) + 0
  npair = t(aa) %*% bb
  P <- 2 * (1 - pt(abs(h) * sqrt(npair - 2)/sqrt(1 - 
                                                   h * h), npair - 2))
  P[abs(h) == 1] <- 0
  P[is.na(P)] <-1
  list(estimate=h, P = P)
}

