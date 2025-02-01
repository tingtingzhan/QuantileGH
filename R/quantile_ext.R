


#' @title Variance-Covariance of Quantiles
#' 
#' @description
#' To compute the variance-covariance matrix of quantiles 
#' based on Theorem 1 and 2 of Mosteller (1946).
#' 
#' @param p \link[base]{numeric} \link[base]{vector}, 
#' cumulative probabilities at the given quantiles
#' 
#' @param d \link[base]{numeric} \link[base]{vector}, 
#' probability densities at the given quantiles
#' 
#' @details 
#' 
#' The end user should make sure no density too close to 0 is included in argument `d`.
#' 
#' Function [quantile_vcov] must not be used in a compute-intensive way.
#' 
#' @returns 
#' 
#' Function [quantile_vcov] returns the variance-covariance \link[base]{matrix} of quantiles.
#' 
#' @references 
#' Frederick Mosteller. On Some Useful "Inefficient" Statistics (1946).
#' \doi{10.1214/aoms/1177730881}
#' 
#' @keywords internal
#' @export
quantile_vcov <- function(p, d) {
  # do the check on d=0 in [QLMDe]
  if (anyNA(p) || anyNA(d)) stop('no NA allowed in probability nor density')
  if ((n <- length(p)) != length(d)) stop('p and d must match in length')
  
  fs <- tcrossprod(d, d) # 'matrix'
  p_c <- array(p, dim = c(n,n)) # 'p on cols'
  p_r <- t.default(p_c) # 'p on rows'
  p_min <- pmin.int(p_r, p_c) # vector!
  p_max <- pmax.int(p_r, p_c)
  vv <- p_min * (1 - p_max) / fs # back to 'matrix'
  return(vv)
}