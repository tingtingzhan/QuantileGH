


#' @title Forward Selection of the Number of Components \eqn{K}
#' 
#' @description
#' 
#' To compare \eqn{gh}-parsimonious models of Tukey \eqn{g}-&-\eqn{h} mixtures with different number of components \eqn{K} 
#' (up to a user-specified \eqn{K_\text{max}})
#' and select the optimal number of components.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, observations
#' 
#' @param distname,data.name \link[base]{character} scalars, 
#' see parameters of the same names in function [QLMDe()]
#' 
#' @param direction \link[base]{character} scalar, direct of selection in function [step_fmx()], 
#' either `'forward'` (default) or `'backward'`
#' 
#' @param test \link[base]{character} scalar, criterion to be used, either 
#' Akaike's information criterion \link[stats]{AIC}, or 
#' Bayesian information criterion \link[stats]{BIC} (default).
#' 
#' @param Kmax \link[base]{integer} scalar \eqn{K_\text{max}}, 
#' maximum number of components to be considered. Default `3L`
#' 
#' @param ... additional parameters
#' 
#' @details 
#' 
#' Function [QLMDe_stepK()] compares the \eqn{gh}-parsimonious models with different number of components \eqn{K},
#' and selects the optimal number of components using BIC (default) or AIC.
#' 
#' The forward selection starts with finding the \eqn{gh}-parsimonious model (via function [step_fmx()])
#' at \eqn{K = 1}.   
#' Let the current number of component be \eqn{K^c}.  
#' We compare the \eqn{gh}-parsimonious models of \eqn{K^c+1} and \eqn{K^c} component, respectively,
#' using BIC or AIC. 
#' If \eqn{K^c} is preferred, then the forward selection is stopped, and \eqn{K^c} is considered the 
#' optimal number of components.
#' If \eqn{K^c+1} is preferred, then
#' the forward selection is stopped if \eqn{K^c+1=K_{max}},
#' otherwise update \eqn{K^c} with \eqn{K_c+1} and repeat the previous steps.
#' 
#' 
#' @returns 
#' 
#' Function [QLMDe_stepK()] returns an object of S3 class `'stepK'`, 
#' which is a \link[base]{list} of selected models (in reversed order) with attribute(s)
#' `'direction'` and
#' `'test'`.
#' 
#' @examples 
#' data(bmi, package = 'mixsmsn')
#' hist(x <- bmi[[1L]])
#' \donttest{QLMDe_stepK(x, distname = 'GH', Kmax = 2L)}
#' 
#' @export
QLMDe_stepK <- function(
    x, distname = c('GH', 'norm'), data.name = deparse1(substitute(x)), 
    Kmax = 3L, 
    test = c('BIC', 'AIC'), 
    direction = c('forward', 'backward'),
    ...
) {
  
  if ('constraint' %in% names(list(...))) stop('do not specify `constraint`')
  test <- match.arg(test)
  direction <- match.arg(direction)
  distname <- match.arg(distname)
  
  K <- 1L
  
  run0 <- QLMDe(x, distname = distname, data.name = data.name, K = K, 
                constraint = switch(direction, forward = all_constraints_(distname = distname, K = K)), ...)
  mods <- step_fmx(run0, test = test, direction = direction)[1L] # list

  while ((K + 1L) <= Kmax) { # `K + 1` vs. `K`
    new0 <- QLMDe(x, distname = distname, data.name = data.name, K = K + 1L, 
                  constraint = switch(direction, forward = all_constraints_(distname = distname, K = K + 1L)), ...)
    new_ <- step_fmx(new0, test = test, direction = direction)[1L] # list
    
    tval <- lapply(c(new_, mods[1L]), FUN = match.fun(test)) # test value
    if (tval[[2]] <= tval[[1L]] + 1e-07) break # smaller K is selected 

    K <- K + 1L
    mods <- c(new_, mods) # newest fit *in the first element* 
  }
  
  attr(mods, which = 'direction') <- direction 
  attr(mods, which = 'test') <- test 
  class(mods) <- 'stepK'
  return(mods)
}


#' @importFrom fmx print.fmx npar.fmx getTeX
#' @export
print.stepK <- function(x, ...) {
  print.fmx(x[[1L]])
  
  test <- attr(x, which = 'test', exact = TRUE)
  tb <- data.frame( # this is *not* an 'anova' table!!
    '# Parameter' = vapply(x, FUN = npar.fmx, FUN.VALUE = 0L), 
    test = vapply(x, FUN = match.fun(test), FUN.VALUE = 0), 
    row.names = vapply(x, FUN = getTeX, print_K = TRUE, FUN.VALUE = ''), 
    check.names = FALSE)
  names(tb)[2L] <- test
  print.data.frame(tb)
  
  cat('\nUse ', deparse1(substitute(x)), '[[1]] to obtain the selected model\n\n', sep = '')
}






