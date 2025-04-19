

#' @title Tukey \eqn{g}-&-\eqn{h} Transformation
#' 
#' @description
#' Function [z2gh()] (not compute intensive) transforms standard normal quantiles \eqn{z_q} to Tukey \eqn{g}-&-\eqn{h} quantiles.
#' 
#' Function [gh2z()] (compute intensive!!) transforms Tukey \eqn{g}-&-\eqn{h} quantiles to standard normal quantiles \eqn{z_q}.
#'  
#' @param z \link[base]{double} scalar or \link[base]{vector}, standard normal quantiles \eqn{z_q}
#' 
#' @param q \link[base]{double} \link[base]{vector}, Tukey \eqn{g}-&-\eqn{h} quantiles \eqn{q}
#' with location parameter \eqn{A=0} and scale parameter \eqn{B=1}
#' 
#' @param g,h \link[base]{double} scalars
#' 
#' @param ... parameters of function [vuniroot2()], other than `interval`
#' 
#' @details
#' Domain of standard normal quantile \eqn{z} to search for function [gh2z()] is hard coded as 
#' \eqn{(-8.3, 8.3)}, because `stopifnot(identical(pnorm(8.3), 1))`.
#' 
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
gh2z <- function(q, g = 0, h = 0, ...) {
  
  g0 <- (g == 0)
  h0 <- (h == 0)
  
  if (g0 && h0) return(q)
  
  out <- q
  qok <- is.finite(q)
  
  int <- c(-8.3, 8.3)
  # identical(pnorm(8.3), 1) |> stopifnot() # even 8.29 does not work!
  
  if (!g0 && h0) { # has bound but also has explicit form!
    egz <- q[qok]*g + 1
    if (any(id <- (egz <= 0))) {
      out[qok][id] <- if (g > 0) int[1L] else int[2L]
    }
    out[qok][!id] <- log(egz[!id]) / g
    return(out)
  }
  
  if (!g0 && !h0) {
    out[qok] <- vuniroot2(y = q[qok]*g, f = \(z) expm1(g*z) * exp(h*z^2/2), interval = int, ...)
    return(out)
  }
  
  
  if (g0 && !h0) {
    out[qok] <- vuniroot2(y = q[qok], f = \(z) z * exp(h*z^2/2), interval = int, ...)
    return(out)
  }
  
}





