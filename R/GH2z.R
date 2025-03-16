

#' @title Inverse of Tukey \eqn{g}-&-\eqn{h} Transformation
#' 
#' @description
#' To transform Tukey \eqn{g}-&-\eqn{h} quantiles to standard normal quantiles.
#' 
#' @param q \link[base]{double} \link[base]{vector}, quantiles \eqn{q}
#' 
#' @param q0 (optional) \link[base]{double} \link[base]{vector}, 
#' standardized quantiles \eqn{q_0=(q-A)/B}
#' 
#' @param A,B (optional) \link[base]{double} *scalars*, location and scale parameters of 
#' Tukey \eqn{g}-&-\eqn{h} transformation.  Ignored if `q0` is provided.
#' 
#' @param ... parameters of internal helper function [.GH2z()]
#' 
#' @details
#' Unfortunately, function [GH2z()], the inverse of Tukey \eqn{g}-&-\eqn{h} transformation, 
#' does not have a closed form and needs to be solved numerically.
#' 
#' For compute intensive jobs, use internal helper function [.GH2z()].
#' 
#' 
#' @returns 
#' Function [GH2z()] returns a \link[base]{double} \link[base]{vector} of the same length as input `q`.
#' 
#' @examples
#' z = rnorm(1e3L)
#' all.equal.numeric(.GH2z(z2GH(z, g = .3, h = .1), g = .3, h = .1), z)
#' all.equal.numeric(.GH2z(z2GH(z, g = 0, h = .1), g = 0, h = .1), z)
#' all.equal.numeric(.GH2z(z2GH(z, g = .2, h = 0), g = .2, h = 0), z)
#' 
#' @keywords internal
#' @export
GH2z <- function(q, q0 = (q - A)/B, A = 0, B = 1, ...) {
  # ?base::is.finite finds finite AND non-missing; as fast as `rep(TRUE, times = nq)` (where nq = length(q))
  if (!length(q0)) return(numeric()) # required by ?fitdistrplus::fitdist
  out <- q0
  qok <- is.finite(q0)
  out[qok] <- .GH2z(q0 = q0[qok], ...)
  return(out)
}
# @note
# Inspired by `OpVaR:::gh_inv` (package archived)





# internal workhorse of [GH2z]
# inverse of Tukey GH transformation; compute intensive!!!
#' @rdname TukeyGH_helper
#' @export
.GH2z <- function(
    q, q0 = (q - A)/B, # `q` and `q0` both finite AND non-missing
    A = 0, B = 1, g = 0, h = 0, # all scalar (`A` and `B` not needed if `q0` is provided)
    interval = c(-15, 15), # smaller `interval` for \code{QLMDe} algorithm (support of standard normal distribution)
    tol = .Machine$double.eps^.25, maxiter = 1000
) {
  
  #if (!length(q0)) return(numeric()) # required by \link[fitdistrplus]{fitdist}
  g0 <- (g == 0)
  h0 <- (h == 0)
  out <- q0
  #interval <- t.default(array(interval, dim = c(2L, length(q0))))
  
  # bound issue only in [dGH], not [qfmx]
  
  if (!g0 && !h0) { # most likely to happen in ?stats::optim; put in as the first option to save time
    out[] <- vuniroot2(y = q0 * g, f = function(z) expm1(g*z) * exp(h * z^2/2), interval = interval, tol = tol, maxiter = maxiter)
    # very small `h` would cause bound-issue
    return(out)
  }
  
  if (g0 && h0) return(q0)
  
  if (!g0 && h0) { # has bound but also has explicit form!
    egz <- q0 * g + 1
    if (any(id <- (egz <= 0))) {
      out[id] <- if (g < 0) Inf else -Inf
    }
    out[!id] <- log(egz[!id]) / g
    return(out)
  }
  
  if (g0 && !h0) { # wont have the bound issue if g0
    out[] <- vuniroot2(y = q0, f = function(z) z * exp(h * z^2/2), interval = interval, tol = tol, maxiter = maxiter)
    return(out)
  }
  
}




