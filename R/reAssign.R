
#' @title Re-Assign Observations Trimmed Prior to Trimmed \eqn{k}-Means Clustering
#' 
#' @description 
#' Re-assign the observations,
#' which are trimmed in the trimmed \eqn{k}-means algorithm,
#' back to the closest cluster as determined by the smallest 
#' Mahalanobis distance. 
#' 
#' @param x a \link[tclust]{tkmeans} object
#' 
#' @param ... potential parameters, currently not in use.
#' 
#' @details 
#' 
#' Given the \link[tclust]{tkmeans} input, 
#' the \link[stats]{mahalanobis} distance is computed between each trimmed observation 
#' and each cluster.
#' Each trimmed observation is assigned to the closest cluster (i.e., with the smallest Mahalanobis distance). 
#'
#' @returns 
#' Function [reAssign.tkmeans] returns an `'reAssign_tkmeans'` object, 
#' which inherits from \link[tclust]{tkmeans} class. 
#' 
#' @note 
#' Either \link[stats]{kmeans} or \link[tclust]{tkmeans} is slow for big `x`.
#' 
#' 
#' @examples 
#' library(tclust)
#' data(geyser2)
#' clus = tkmeans(geyser2, k = 3L, alpha = .03)
#' plot(clus, main = 'Before Re-Assigning')
#' plot(reAssign(clus), main = 'After Re-Assigning')
#' 
#' @name reAssign
#' @export
reAssign <- function(x, ...) UseMethod('reAssign')

#' @rdname reAssign
#' @importFrom stats cov mad
#' @export reAssign.tkmeans
#' @export
reAssign.tkmeans <- function(x, ...) {
  obs <- x$par$x
  n <- x$int$dim[1L]
  
  cluster <- as.integer(x$cluster) # was 'numeric'
  if (length(cluster) != n) stop('tkmeans: I dont understand?')
  tid <- !cluster # trimmed id
  cluster[tid] <- NA_integer_
  k_seq <- seq_len(x$k)
  k_seq_c <- as.character.default(k_seq)

  fclus <- cluster
  attr(fclus, 'levels') <- k_seq_c
  class(fclus) <- 'factor'
  
  # distance from each of the k clusters
  if ((d <- x$int$dim[2L]) == 1L) {
    obss <- split.default(obs, f = fclus)
    centers <- c(x$centers)
    x_mad <- vapply(k_seq, FUN = function(i) mad(obss[[i]], center = centers[i]), FUN.VALUE = 0)
    tobs <- obs[tid, , drop = TRUE]
    tdist <- t.default(abs(tcrossprod(1/x_mad, tobs) - centers/x_mad)) # not compute intensive
  } else {
    clus_id <- split.default(seq_len(n), f = fclus)
    tobs_t <- t.default(obs[tid, , drop = FALSE]) # trimmed obs
    t_seq <- seq_len(dim(tobs_t)[2L])
    tdist <- do.call(cbind, args = lapply(k_seq, FUN = function(i) {
      ivv <- cov(obs[clus_id[[i]], , drop = FALSE])
      invcov <- chol2inv(chol.default(ivv))
      tobs_i <- tobs_t - x$centers[,i]
      (t.default(tobs_i) %*% invcov %*% tobs_i)[cbind(t_seq, t_seq)]
    }))
  } 
  cluster[tid] <- max.col(-tdist)
    
  out <- x
  out$cluster <- cluster # needs 'integer' or 'double' (not 'factor') for ?tclust:::plot.tkmeans
  out$size <- tabulate(cluster)
  out$weights <- out$size / n
  out$par$alpha <- 0
  # additional elements needs to be changed?
  attr(out, 'class') <- c('reAssign_tkmeans', 'tkmeans')
  return(out)
}