

#' @title Quantile Least Mahalanobis Distance estimates
#' 
#' @description 
#' 
#' The quantile least Mahalanobis distance algorithm estimates the parameters of 
#' single-component or finite mixture distributions   
#' by minimizing the Mahalanobis distance between the vectors of sample and theoretical quantiles.
#' See [QLMDp] for the default selection of probabilities at which the sample and theoretical quantiles are compared.
#' 
#' The default initial values are estimated based on trimmed \eqn{k}-means 
#' clustering with re-assignment.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, the one-dimensional observations.
#' 
#' @param data.name \link[base]{character} scalar, name for the observations for user-friendly print out.
#' 
#' @param distname \link[base]{character} scalar, name of mixture distribution to be fitted.  Currently supports `'norm'` and `'GH'`.
#' 
#' @param K \link[base]{integer} scalar, number of components (e.g., must use `2L` instead of `2`).
#' 
#' @param probs \link[base]{numeric} \link[base]{vector}, percentiles at where the sample and theoretical quantiles are to be matched.
#' See function [QLMDp] for details.
#' 
#' @param init \link[base]{character} scalar for the method of initial values selection, 
#' or an \linkS4class{fmx} object of the initial values. 
#' See function [fmx_hybrid] for more details.
#' 
#' @param constraint \link[base]{character} \link[base]{vector}, parameters (\eqn{g} and/or \eqn{h} for Tukey \eqn{g}-&-\eqn{h} mixture) to be set at 0.  
#' See function [fmx_constraint] for details.
#' 
#' @param tol,maxiter see function \link[TukeyGH77]{vuniroot2}
#' 
#' @param ... additional parameters of \link[stats]{optim}
#' 
#' @details
#' 
#' Quantile Least Mahalanobis Distance estimator fits a single-component or finite mixture distribution 
#' by minimizing the Mahalanobis distance between
#' the theoretical and observed quantiles,
#' using the empirical quantile variance-covariance matrix [quantile_vcov].
#' 
#' @returns 
#' 
#' Function [QLMDe] returns an \linkS4class{fmx} object.
#' 
#' 
#' 
#' @examples 
#' data(bmi, package = 'mixsmsn')
#' hist(x <- bmi[[1L]])
#' \donttest{QLMDe(x, distname = 'GH', K = 2L)}
#' 
#' @seealso [fmx_hybrid]
#' @importFrom fmx distArgs
#' @importFrom fmx Kolmogorov_fmx CramerVonMises_fmx KullbackLeibler_fmx
#' @importFrom fmx approxdens fmx2dbl dbl2fmx pmlogis_first dfmx qfmx fmx_constraint user_constraint
#' @importFrom stats ecdf optim pnorm qnorm quantile
#' @importFrom TukeyGH77 vuniroot2 qGH .GH2z
#' @export
QLMDe <- function(
  x, distname = c('GH', 'norm', 'sn'), K, data.name = deparse1(substitute(x)),
  constraint = character(),
  probs = QLMDp(x = x),
  init = c('logLik', 'letterValue', 'normix'),
  tol = .Machine$double.eps^.25, maxiter = 1000,
  ...
) {
  
  distname <- match.arg(distname)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || K <= 0L) stop('number of component must be length-1 positive integer')
  if (!is.integer(K)) stop('number of component must be length-1 positive integer (e.g., use integer `2L` instead of numeric/double `2` for 2-component mixture model)')
  
  if (!is.vector(x, mode = 'double')) stop('x must be double vector')
  if (anyNA(x)) stop('do not allow NA_real_ in observations `x`')
  if (!is.character(data.name) || length(data.name) != 1L || anyNA(data.name) || !nzchar(data.name)) stop('data.name must be length-1 character')
  
  interval <- QLMDe_interval(x = x, ...)
  
  if (anyNA(probs)) stop('do not allow NA_real_ in `probs`')
  probs <- sort.int(unique_allequal(probs))
  
  if (missing(init)) {
    init <- switch(match.arg(init), logLik = {
      fmx_hybrid(x, test = 'logLik', distname = distname, K = K, constraint = constraint)
    }, letterValue = {
      fmx_cluster(x, distname = distname, K = K, constraint = constraint)
    }, normix = {
      fmx_normix(x, distname = distname, K = K)
    })
  }
  if (!inherits(init, what = 'fmx')) stop('`init` must be \'fmx\' object')
  if ((init@distname != distname) || (dim(init@pars)[1L] != K)) stop('`init` is not a ', distname, ' ', K, '-component fit.')
  
  q_init <- qfmx(p = probs, interval = interval, distname = distname, K = K, pars = init@pars, w = init@w)
  if (any(id1 <- is.infinite(q_init))) {
    if (all(id1)) stop('starting values too far away?')
    probs <- probs[!id1]
  }
  q_obs <- quantile(x, probs = probs) # observed quantiles, constant in ?stats::optim
  
  x_epdf <- approxdens(x)
  # Tingting is not sure whether ?stats::approx *and* ?stats::approxfun will make `x_epdf(q_obs) = 0` more likely
  d_obs <- x_epdf(q_obs) # observed density evaluated at `q_obs`
  if (anyNA(d_obs)) stop('do not allow NA_real_ empirical density')
  tol <- sqrt(sqrt(.Machine$double.eps))
  if (all(d0 <- (abs(d_obs) < tol))) stop('must have at least one positive density') # `d_obs` should always be positive, but ?base::abs should not hurt
  if (any(d0)) {
    probs <- probs[!d0] # quantiles, where empirical densities are too close to 0, are excluded from the calculation 
    d_obs <- d_obs[!d0]
    q_obs <- q_obs[!d0]
  }
  
  npar <- K * switch(distname, norm = 2L, GH = 4L) + (K - 1L)
  if (length(probs) < npar) {
    stop('Using ', length(probs), ' matching-quantiles to estimate a mixture distribution with ', npar, ' independent parameters. ', 
         'Try increasing `probs` (see ?QLMDp for detail).')
  }
  
  qvv <- quantile_vcov(p = probs, d = d_obs) 
  qvv_inv <- chol2inv(chol.default(qvv))
  
  parRun <- fmx2dbl(init)
  id_constr <- user_constraint(constraint, distname = distname, K = K)
  # the constraint is from user-specification, instead of from Hoaglin's starting value (i.e. h=0 for negative slope)
  has_constr <- (length(id_constr) > 0L)
  
  # ?stats::optim does not allow Inf in `par`
  par_init <- if (has_constr) parRun[-id_constr] else parRun
  if (any(is.infinite(par_init))) {
    par_init[(par_init < 0) & is.infinite(par_init)] <- -5 
    # .Machine$double.min.exp too small (exp(h) trapped at 0) # exp(-5) close to 0 enough
    # e.g., if h=0 was given in Hoaglin's starting value (due to negative slope), but not as a user-specified constraint,
    # I will use exp(-5)=0.0067 as the starting value, which is passed into ?stats::optim
    if (any((par_init > 0) & is.infinite(par_init))) stop('+Inf parameter indicates ???')
  }
  
  # or use this: https://stackoverflow.com/questions/52552143/how-to-save-the-coefficients-for-each-optim-iteration
  max_return <- .Machine$double.xmax
  
  fn <- if (K == 1L) {
    
    switch(distname, norm = function(x) {
      q <- qnorm(probs, mean = x[1L], sd = exp(x[2L]), lower.tail = TRUE, log.p = FALSE)
      if (any(is.infinite(q))) stop('1-comp normal, infinite `q` should not be returned from any `probs` between 0 to 1, see `qnorm(.Machine$double.eps)`') # return(max_return)
      return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      
    }, GH = function(x) {
      if (has_constr) parRun[-id_constr] <- x else parRun <- x
      q <- qGH(probs, A = parRun[1L], B = exp(parRun[2L]), g = parRun[3L], h = exp(parRun[4L]), lower.tail = TRUE, log.p = FALSE)
      if (any(is.infinite(q))) return(max_return) # stop('interval problem again?') #
      return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      
    })
    
  } else {
    
    Kseq <- seq_len(K)
    Kseq1 <- seq_len(K - 1L)

    switch(distname, norm = {
      id_w <- 2L*K + Kseq1
      function(x) {
        .pM <- array(x[seq_len(2L*K)], dim = c(K, 2L)) # only first 2K elements
        t_w <- t.default(pmlogis_first(x[id_w]))
        sdinv <- 1 / exp(.pM[,2L])
        eff <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L]))) * sdinv
        q <- vuniroot2(y = probs, f = function(q) { # essentially [pfmx]
          z <- tcrossprod(sdinv, q) - eff
          c(t_w %*% pnorm(z))
        }, interval = interval)
        if (any(is.infinite(q))) return(max_return) # stop('interval problem again?') # 
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }
      
    }, GH = {
      id_w <- 4L*K + Kseq1
      function(x) {
        if (has_constr) parRun[-id_constr] <- x else parRun <- x
        .pM <- array(parRun[seq_len(4L*K)], dim = c(K, 4L)) # only first 4K elements
        t_w <- t.default(pmlogis_first(parRun[id_w]))
        g <- .pM[,3L]
        h <- exp(.pM[,4L])
        sdinv <- 1 / exp(.pM[,2L])
        eff <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L]))) * sdinv
        q <- vuniroot2(y = probs, f = function(q) { # essentially [pfmx]
          z <- q0 <- tcrossprod(sdinv, q) - eff
          for (i in Kseq) z[i,] <- .GH2z(q0 = q0[i,], g = g[i], h = h[i], tol = tol, maxiter = maxiter)
          c(t_w %*% pnorm(z))
        }, interval = interval, tol = tol, maxiter = maxiter)
        if (any(is.infinite(q))) return(max_return)
        #if (any(is.infinite(q))) {
        #  print(.pM)
        #  print(probs)
        #  print(q)
        #  print(interval)
        #  stop('interval problem again?')
        #}
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }
      
    })
  }
  
  #repeat {
    y <- optim(par = par_init, fn = fn, ...)
    if (FALSE) {
      fn(par_init)
      fn(y$par)
    }
    if (isTRUE(all.equal.numeric(y$par, par_init))) {
      #print(init)
      # ret <<- init; x <<- x;
      stop('stats::optim not working (most likely due to poor starting value)')
      # autoplot.fmx(ret, data = x)
      # QLMDe(x, distname = 'GH', K = 2L, init = ret)
    }
    #cat('run stats::optim again..\n')
    # sometimes ?stats::optim will not move.  Weird
  #}
  
  
  if (has_constr) {
    parRun[-id_constr] <- y$par
  } else parRun <- y$par
  tmp <- dbl2fmx(x = parRun, K = K, distname = distname)
  colnames(tmp$pars) <- distArgs(distname = distname) # will remove after \pkg{fmx} update

  # variance-covariance of internal estimates
  #n <- length(x)
  # `q`: theoretical quantiles
  q <- qfmx(probs, distname = distname, K = K, pars = tmp$pars, w = tmp$w, interval = interval, lower.tail = TRUE, log.p = FALSE)
  # `d`: densities at theoretical (\hat{\theta}) quantiles
  d <- dfmx(q, distname = distname, K = K, pars = tmp$pars, w = tmp$w, log = FALSE)
  # {density at \hat{\theta}} and {kernel density} may have different {=0} status.
  tol <- sqrt(sqrt(.Machine$double.eps))
  if (all(d0 <- (abs(d) < tol))) {
    #cat('malformed estimates with all-0 densities\n')
    #return(invisible())
    stop('do not allow this to happen')
  }
  if (any(d0)) {
    p1 <- probs[!d0]
    q1 <- q[!d0]
    d1 <- d[!d0]
  } else {
    p1 <- probs
    q1 <- q
    d1 <- d
  }
  .meat <- quantile_vcov(p = p1, d = d1) # V_(\hat{theta})
  # in theary, we should use V_{true theta}, but no one knows true theta in practice
  # so I am using V_{\hat_theta}
  # now I want to use V_{empirical} (`@quantile_vv` is no longer a slot of \linkS4class{fmx}), can I?
  q_gr <- qfmx_gr(probs = p1, distname = distname, K = K, pars = tmp$pars, w = tmp$w)
  if (!length(q_gr)) {
    int_vv <- array(NA_real_, dim = c(0L, 0L)) # exception handling
  } else {
    #.bread <- tryCatch(crossprod_inv(q_gr) %*% t.default(q_gr), error = as.null.default)
    .bread <- tryCatch(solve.default(crossprod(q_gr)) %*% t.default(q_gr), error = as.null.default)
    if (!length(.bread)) {
      int_vv <- array(NA_real_, dim = c(0L, 0L)) # exception handling
    } else {
      int_vv <- .bread %*% tcrossprod(.meat, .bread) / length(x) # interval variance-covariance
      if (anyNA(int_vv)) stop('do not allow NA in `int_vv`')
      if (any(diag(int_vv) < 0)) stop('diagonal terms of VarCov < 0 ?!')
      dimnames(int_vv) <- list(dimnames(q_gr)[[2L]], dimnames(q_gr)[[2L]])
    }
  }
  # end of vv
  
  
  ret <- new(
    Class = 'fmx',
    data = x, data.name = data.name,
    distname = distname, pars = tmp$pars, w = tmp$w,
    #quantile_vv = qvv,
    vcov_internal = int_vv
    #probs = probs#,
    # init = init,
    # optim = y
  )
  
  ret@Kolmogorov <- Kolmogorov_fmx(object = ret)
  ret@CramerVonMises <- CramerVonMises_fmx(object = ret)
  ret@KullbackLeibler <- KullbackLeibler_fmx(object = ret)
  
  attr(ret, which = 'probs') <- probs # needed in [drop1.fmx] and [add1.fmx]
  
  if (!setequal(attr(fmx_constraint(ret), which = 'user', exact = TRUE), constraint)) {
    #message('not handling constrants correctly')
    # indicates ?stats::optim did not move
    #return(invisible())
    stop('should not happen')
  }
  
  return(ret)
  
}



QLMDe_interval <- function(x, extend_interval = 10, ...) {
  xmin <- min(x)
  xmax <- max(x)
  xdiff <- xmax - xmin
  if (length(extend_interval) != 1) stop('`extend_interval` must be length-1')
  # interval <- c(min(x), max(x)) # old
  c(xmin - extend_interval * xdiff, xmax + extend_interval * xdiff)
}



# gradient of ?qfmx with respect to the *unconstraint* parameters.
#' @importFrom fmx fmx2dbl fmx_constraint qfmx_interval
#' @importFrom stats setNames numericDeriv
qfmx_gr <- function(
  dist, # can be missing
  probs = stop('must provide `probs`'), 
  # constant given `dist`
  distname = dist@distname, pars = dist@pars, K = dim(pars)[1L], w = dist@w,
  interval = qfmx_interval(distname = distname, pars = pars, K = K, w = w, probs = c(1e-5, 1-1e-5)),
  ...
) {
  x_skeleton <- fmx2dbl(distname = distname, pars = pars, K = K, w = w)
  has_constr <- (length(id_constr <- fmx_constraint(distname = distname, K = K, pars = pars)) > 0L)
  # ?fmx2dbl may have log(d) being -Inf; see ?fmx2dbl
  x_dbl <- if (has_constr) x_skeleton[-id_constr] else x_skeleton
  x_nm <- names(x_dbl) # stopifnot(all(make.names(x_nm) == x_nm)); # unicode symbols considered legal :)
  if (any(is.infinite(x_dbl))) {
    # stop('wont happen since implementation of constraint') # actually could still happen
    return(invisible())
  }

  if (!is.numeric(interval) || length(interval) != 2L || !all(is.finite(interval))) stop('`interval` must not contain NA nor Inf')

  x_sbl <- lapply(x_nm, FUN = as.symbol)
  qfmx_fn <- as.function.default(c(setNames(rep(x = alist(. = ), times = length(x_nm)), nm = x_nm), as.call(c(
    quote(`{`), 
    quote(x <- x_skeleton),
    if (has_constr) {
      as.call(c(quote(`<-`), quote(x[-id_constr]), as.call(c(quote(c), x_sbl))))
    } else as.call(c(quote(`<-`), quote(x), as.call(c(quote(c), x_sbl)))),
    quote(fx <- dbl2fmx(x = x, distname = distname, K = K, argnm = NULL)),
    quote(qfmx(p = probs, pars = fx$pars, w = fx$w, interval = interval, distname = distname, K = K))
  ))))
  
  ret <- tryCatch(expr = {
    # constants `x_skeleton`, `probs`, `interval`, `distname`, `K` does not need to be assigned to `rho` :))
    attr(numericDeriv(
      expr = as.call(c(quote(qfmx_fn), x_sbl)), 
      theta = x_nm, 
      central = TRUE, # needed, otherwise could (very rarely) have gradient for `logitw2` all 0
      rho = list2env(as.list.default(x_dbl), parent = environment())), which = 'gradient', exact = TRUE)
  }, error = function(e) { # very rare
    cat('stats::numericDeriv in `qfmx_gr` error.\n')
    #array(NA_real_, dim = c(length(probs), length(x_nm)))
    return(invisible())
  })
  if (length(ret)) dimnames(ret)[[2L]] <- x_nm
  return(ret)
}




