
##########################################
# Modify tzhInternal/allequal.R
##########################################


#' @title Determine Nearly-Equal Elements
#' 
#' @description 
#' Determine nearly-equal elements and extract non-nearly-equal elements in a \link[base]{double} \link[base]{vector}.
#' 
#' @param x \link[base]{double} \link[base]{vector}
#' 
#' @param ... additional parameters of function [allequal_o_()]
#' 
#' @returns 
#' 
#' Function [duplicated_allequal()] returns a \link[base]{logical} \link[base]{vector} of the same length as the input vector,
#' indicating whether each element is nearly-equal to any of the previous elements.  
#' 
#' Function [unique_allequal()] returns the non-nearly-equal elements in the input vector.
#' 
#' @examples 
#' (x = c(.3, 1-.7, 0, .Machine$double.eps))
#' duplicated.default(x) # not desired
#' unique.default(x) # not desired
#' duplicated_allequal(x)
#' unique_allequal(x)
#' unique_allequal(x, tol = .Machine$double.eps/2)
#' @seealso \link[base]{duplicated.default} \link[base]{unique.default}
#' @name allequal
#' @keywords internal
#' @export
unique_allequal <- function(x, ...) x[!duplicated_allequal(x, ...)]


#' @rdname allequal
#' @export
duplicated_allequal <- function(x, ...) {
  x <- unclass(x)
  if (!is.vector(x, mode = 'double')) stop('input must be double vector')
  nx <- length(x)
  id <- allequal_o_(target = x, current = x, ...)
  id[lower.tri(id, diag = TRUE)] <- FALSE # super smart!
  return(.colSums(id, m = nx, n = nx, na.rm = FALSE) > 0) # takes care of 1st element too 
}





#' @title Test if Two \link[base]{double} Vectors are Element-Wise (Nearly) Equal
#' 
#' @description
#' 
#' Test if two \link[base]{double} \link[base]{vector}s are element-wise (nearly) equal.
#' 
#' @param current length-\eqn{n_c} \link[base]{double} \link[base]{vector}, 
#' the value(s) to be compared with `target`, missing value not allowed
#' 
#' @param target length-\eqn{n_t} \link[base]{double} \link[base]{vector}, 
#' the target value(s), missing value not allowed
#' 
#' @param tolerance positive \link[base]{double} scalar, default `sqrt(.Machine$double.eps)` 
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details 
#' 
#' Function [allequal_o_()] is different from \link[base]{all.equal.numeric}, such that 
#' \itemize{
#' \item{only compares between \link[base]{double}, not \link[base]{complex}, values}
#' \item{element-wise comparison is performed, in a fashion similar to function \link[base]{outer}}
#' \item{a \link[base]{logical} scalar is always returned for each element-wise comparison.}
#' }
#' 
#' @returns 
#' Function [allequal_o_()] returns an \eqn{n_t\times n_c} \link[base]{logical} \link[base]{matrix}
#' indicating whether the length-\eqn{n_c} \link[base]{vector} `current` 
#' is element-wise near-equal to 
#' the length-\eqn{n_t} \link[base]{vector} `target` 
#' within the pre-specified `tolerance`.  
#' 
#' @examples 
#' allequal_o_(target = c(.3, 0), current = c(.3, 1-.7, 0, .Machine$double.eps))
#' @keywords internal
#' @export
allequal_o_ <- function(target, current, tolerance = sqrt(.Machine$double.eps), ...) {
  if (!(nt <- length(target)) || !is.double(target)) stop('illegal `target`')
  if (!(nc <- length(current)) || !is.double(current)) stop('illegal `current`')
  if (anyNA(target) || anyNA(current)) stop('Do not allow missingness in `target` or `current`')
  
  t. <- rep(target, times = nc) # `target` on row
  c. <- rep(current, each = nt) # `current` on col
  # all.equal.numeric(target = t., current = c., ...) # nah, cannot do this..
  
  # object names inspired by ?base::all.equal.numeric
  xy <- abs(t. - c.) # 'absolute'
  xn <- abs(t.)
  if (all(is.finite(xn)) && all(xn > tolerance)) xy <- xy/xn # 'relative'
  ret <- (xy <= tolerance)
  dim(ret) <- c(nt, nc) # put to last step
  return(ret)
}
