

#' @title Quantile Least Mahalanobis Distance Estimator for Tukey \eqn{g}-&-\eqn{h} Mixture
#'
#' @description
#'
#' Tools for simulating and fitting finite mixtures of the 4-parameter Tukey \eqn{g}-&-\eqn{h} distributions. 
#' Tukey \eqn{g}-&-\eqn{h} mixture is highly flexible to model multimodal distributions with variable degree of skewness and kurtosis in the components. 
#' The Quantile Least Mahalanobis Distance estimator [QLMDe] is used for estimating parameters of the finite Tukey \eqn{g}-&-\eqn{h} mixtures.
#' [QLMDe] is an indirect estimator that minimizes the Mahalanobis distance between the sample and model-based quantiles.
#' A backward-forward stepwise model selection algorithm is provided to find
#' \itemize{
#' \item {a parsimonious Tukey \eqn{g}-&-\eqn{h} mixture model, conditional on a given number-of-components; and}
#' \item {the optimal number of components within the user-specified range.}
#' }
#'
# @references none yet
#'
#'
#' @examples
#' # see ?QLMDe
#'
'_PACKAGE'








