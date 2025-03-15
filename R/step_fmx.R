
#' @title Forward Selection of \eqn{gh}-parsimonious Model with Fixed Number of Components \eqn{K}
#' 
#' @description 
#' 
#' To select the \eqn{gh}-parsimonious mixture model, 
#' i.e., with some \eqn{g} and/or \eqn{h} parameters equal to zero,
#' conditionally on a fixed number of components \eqn{K}.
#' 
#' @param object \link[fmx:fmx-class]{fmx} object
#' 
#' @param test \link[base]{character} scalar, criterion to be used, either 
#' Akaike's information criterion \link[stats]{AIC}-like, or 
#' Bayesian information criterion \link[stats]{BIC}-like (default).
#' 
#' @param direction \link[base]{character} scalar, `'forward'` (default) or
#' `'backward'`
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' The algorithm starts with quantile least Mahalanobis distance estimates
#' of either the full mixture of Tukey \eqn{g}-&-\eqn{h} distributions model, or
#' a constrained model (i.e., some \eqn{g} and/or \eqn{h} parameters equal to zero according to the user input).
#' Next, each of the non-zero \eqn{g} and/or \eqn{h} parameters is tested using the likelihood ratio test.
#' If all tested \eqn{g} and/or \eqn{h} parameters are significantly different from zero at the level 0.05
#' the algorithm is stopped and the initial model is considered \eqn{gh}-parsimonious.
#' Otherwise, the \eqn{g} or \eqn{h} parameter with the largest p-value is constrained to zero 
#' for the next iteration of the algorithm.
#' 
#' The algorithm iterates until only significantly-different-from-zero \eqn{g} and \eqn{h} parameters 
#' are retained, which corresponds to \eqn{gh}-parsimonious Tukey \eqn{g}-&-\eqn{h} mixture model.
#' 
#' @returns 
#' 
#' Function [step_fmx()] returns an object of S3 class `'step_fmx'`, 
#' which is a \link[base]{list} of selected models (in reversed order) with attribute(s)
#' `'direction'` and
#' `'test'`.
#' 
#' @seealso 
#' \link[stats]{step}
#' 
#' @importClassesFrom fmx fmx
#' @export
step_fmx <- function(object, test = c('BIC', 'AIC'), direction = c('forward', 'backward'), ...) {
  if (!length(object@data)) return(invisible())
  test <- match.arg(test)
  direction <- match.arg(direction)
  K <- dim(object@pars)[1L]
  obj_start <- QLMDe(x = object@data, distname = object@distname, data.name = object@data.name, K = K, 
                     constraint = switch(direction, forward = all_constraints_(distname = object@distname, K = K)), ...)
  mods <- list(obj_start)
  message('Finding parsimonious mixture at K = ', K)
  repeat {
    tmp <- c(list(mods[[1L]]), # running model as the 1st element
             switch(direction, backward = drop1.fmx, forward = add1.fmx)(object = mods[[1L]], ...))
    o1 <- order(vapply(tmp, FUN = match.fun(test), FUN.VALUE = 0, USE.NAMES = FALSE), decreasing = FALSE)[1L]
    if (o1 == 1L) break # running model is the best
    mods <- c(list(tmp[[o1]]), mods) # new selection appended to 1st index
  }
  attr(mods, which = 'direction') <- direction # I dont think I used this anywhere..
  attr(mods, which = 'test') <- test
  class(mods) <- 'step_fmx'
  return(mods)
}


#' @importFrom fmx print.fmx npar.fmx getTeX
#' @export
print.step_fmx <- function(x, ...) {
  print.fmx(x[[1L]])
  
  test <- attr(x, which = 'test', exact = TRUE)
  tb <- data.frame( # this is *not* an 'anova' table!!
    '# Parameter' = vapply(x, FUN = npar.fmx, FUN.VALUE = 0L), 
    test = vapply(x, FUN = match.fun(test), FUN.VALUE = 0), 
    row.names = vapply(x, FUN = getTeX, FUN.VALUE = ''), 
    check.names = FALSE)
  names(tb)[2L] <- test
  print.data.frame(tb)
  
  cat('\nUse ', deparse1(substitute(x)), '[[1]] to obtain the selected model\n\n', sep = '')
}




all_constraints_ <- function(distname, K) {
  switch(distname, GH = {
    c(outer(c('g', 'h'), seq_len(K), FUN = paste0))
  }, character())
}

