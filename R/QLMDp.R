
#' @title Percentages for Quantile Least Mahalanobis Distance estimation
#' 
#' @description 
#' 
#' A vector of probabilities to be used in Quantile Least Mahalanobis Distance estimation ([QLMDe]). 
#' 
#' @param from,to \link[base]{numeric} scalar, 
#' minimum and maximum of the equidistant (in probability or quantile) probabilities.  
#' Default `.05` and `.95`, respectively
#' 
#' @param length.out non-negative \link[base]{integer} scalar, 
#' the number of the equidistant (in probability or quantile) probabilities. 
#' 
#' @param equidistant \link[base]{character} scalar.
#' If `'prob'` (default), then the probabilities are equidistant.  
#' If `'quantile'`, then the quantiles (of the observations `x`) corresponding to the probabilities are equidistant.
#' 
#' @param extra \link[base]{numeric} \link[base]{vector} of *additional* probabilities, 
#' default `c(.005, .01, .02, .03, .97, .98, .99, .995)`.
#' 
#' @param x \link[base]{numeric} \link[base]{vector} of observations, only used when `equidistant = 'quantile'`.
#' 
#' @details
#' 
#' The default arguments of function [QLMDp] returns the probabilities of 
#' `c(.005, .01, .02, .03, seq.int(.05, .95, length.out = 15L), .97, .98, .99, .995)`.
#' 
#' @returns 
#' 
#' A \link[base]{numeric} \link[base]{vector} of probabilities to be supplied to parameter `p` of 
#' Quantile Least Mahalanobis Distance [QLMDe] estimation).
#' In practice, the length of this probability \link[base]{vector} `p` 
#' must be equal or larger than the number of parameters in the distribution model to be estimated.
#' 
#' @examples 
#' library(fmx)
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(100); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' 
#' # equidistant in probabilities
#' (p1 = QLMDp()) 
#' 
#' # equidistant in quantiles
#' (p2 = QLMDp(equidistant = 'quantile', x = x2)) 
#' @export
QLMDp <- function(
  from = .05, to = .95, length.out = 15L, 
  equidistant = c('prob', 'quantile'),
  extra = c(.005, .01, .02, .03, .97, .98, .99, .995),
  x
) {
  
  if (!is.double(from) || length(from) != 1L || is.na(from) || from <= 0) stop('`from` must be numeric >0')
  if (!is.double(to) || length(to) != 1L || is.na(to) || to >= 1) stop('`to` must be numeric <1')
  if (!is.integer(length.out) || length(length.out) != 1L || anyNA(length.out) || length.out < 0L) stop('N must be len-1 non-negative integer')
  if (!is.double(extra) || anyNA(extra) || any(extra <= 0, extra >= 1)) stop('`extra` must be numeric between (0, 1)') # complatible with len-0

  p <- if (!length.out) extra else {
    c(extra, switch(match.arg(equidistant), prob = {
      #cat('Equidistant in probability; ')
      seq.int(from = from, to = to, length.out = length.out)
    }, quantile = {
      #cat('Equidistant in quantile; ')
      if (missing(x)) stop('Must have `x`')
      qlim <- quantile(x, probs = c(from, to)) # not compute intensive
      ecdf(x)(seq.int(from = qlim[1L], to = qlim[2L], length.out = length.out))
    }))
  }
  if (!length(p)) stop('len-0 p, do not allow')
  sort.int(unique_allequal(p))
}




