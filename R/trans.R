

#' @title Tukey \eqn{g}-&-\eqn{h} Transformation
#' 
#' @description
#' To transform 
#' standard normal quantiles \eqn{z_q} 
#' to 
#' Tukey \eqn{g}-&-\eqn{h} quantiles, 
#' or vise versa.
#'  
#' @param z \link[base]{double} scalar or \link[base]{vector}, standard normal quantiles \eqn{z_q}
#' 
#' @param q \link[base]{double} \link[base]{vector}, Tukey \eqn{g}-&-\eqn{h} quantiles \eqn{q}
#' with location parameter \eqn{A=0} and scale parameter \eqn{B=1}
#' 
#' @param g,h \link[base]{double} scalars
#' 
#' @param interval \link[base]{numeric} \link[base]{length}-`2L` \link[base]{vector} for function [gh2z()],
#' domain of standard normal quantile \eqn{z} to search,
#' default `c(-8.3, 8.3)` as `stopifnot(identical(pnorm(8.3), 1))`.
#' See more from function [vuniroot2()].
#' 
#' @param ... other parameters of function [vuniroot2()]
#' 
#' @details
#' Function [z2gh()] (not compute intensive) transforms standard normal quantiles \eqn{z_q} to Tukey \eqn{g}-&-\eqn{h} quantiles.
#' 
#' Function [gh2z()] (compute intensive!!) transforms Tukey \eqn{g}-&-\eqn{h} quantiles to standard normal quantiles \eqn{z_q}.
#' 
#' @returns
#' Functions [z2gh()] and [gh2z()] both return 
#' \link[base]{double} scalar, \link[base]{vector} or \link[base]{matrix}.
#' 
#' @keywords internal
#' @name tukey_transform
#' @export
z2gh <- function(z, g = 0, h = 0) {
  # `z` must be numeric vector (i.e. not 'matrix'), all arguments `g` and `h` are scalars
  g0 <- (g == 0) # scalar
  h0 <- (h == 0)
  if (g0 && h0) return(z)
  if (g0 && !h0) return(z * exp(h*z^2/2))
  if (!g0 && h0) return(expm1(g*z) / g)
  return(expm1(g*z) / g * exp(h*z^2/2))
}

#' @rdname tukey_transform
#' @export
gh2z <- function(
    q,
    g = 0, h = 0,
    interval = c(-8.3, 8.3),
    ...
) {
  
  g0 <- (g == 0)
  h0 <- (h == 0)
  
  if (g0 && h0) return(q)
  
  out <- q
  qok <- is.finite(q)
  
  if (!g0 && h0) { # has bound but also has explicit form!
    egz <- q[qok]*g + 1
    if (any(id <- (egz <= 0))) {
      out[qok][id] <- if (g > 0) interval[1L] else interval[2L]
    }
    out[qok][!id] <- log(egz[!id]) / g
    return(out)
  }
  
  if (!g0 && !h0) { # most likely to happen in ?stats::optim; put in as the first option to save time
    out[qok] <- vuniroot2(y = q[qok]*g, f = \(z) expm1(g*z) * exp(h*z^2/2), interval = interval, ...)
    # very small `h` would cause bound-issue
    return(out)
  }
  
  
  if (g0 && !h0) { # wont have the bound issue if g0
    out[qok] <- vuniroot2(y = q[qok], f = \(z) z * exp(h*z^2/2), interval = interval, ...)
    return(out)
  }
  
}





