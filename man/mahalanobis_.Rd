% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mahalanobis_.R
\name{mahalanobis_}
\alias{mahalanobis_}
\title{Simpler and Faster Mahalanobis Distance}
\usage{
mahalanobis_(x, center, invcov)
}
\arguments{
\item{x}{\link[base]{numeric} \link[base]{vector}}

\item{center}{\link[base]{numeric} \link[base]{vector}, mean \eqn{\mathbf{\mu}}}

\item{invcov}{\link[base]{numeric} \link[base]{matrix}, \emph{inverted} variance-covariance \eqn{\mathbf{\Sigma}}}
}
\value{
Function \code{\link[=mahalanobis_]{mahalanobis_()}} returns a \link[base]{numeric} scalar.
}
\description{
Simpler and faster \link[stats]{mahalanobis} distance.
}
\keyword{internal}
