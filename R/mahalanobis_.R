

#' @title Simpler and Faster Mahalanobis Distance
#' 
#' @description 
#' 
#' Simpler and faster \link[stats]{mahalanobis} distance.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}
#' 
#' @param center \link[base]{numeric} \link[base]{vector}, mean \eqn{\mathbf{\mu}}
#' 
#' @param invcov \link[base]{numeric} \link[base]{matrix}, *inverted* variance-covariance \eqn{\mathbf{\Sigma}}
#' 
#' @returns 
#' 
#' Function [mahalanobis_()] returns a \link[base]{numeric} scalar.
#' 
#' @keywords internal
#' @export
mahalanobis_ <- function(x, center, invcov) {
  # if (!is.vector(x, mode = 'double')) stop('x must be double vector') # speed
  # if (!is.vector(center, mode = 'double')) stop('center must be double vector') # speed
  x0 <- x - center
  c(crossprod(x0, invcov) %*% x0)
}


