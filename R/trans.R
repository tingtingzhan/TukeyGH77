

#' @title Tukey \eqn{g}-&-\eqn{h} Transformation
#' 
#' @description
#' Function [z2gh()] (not compute intensive) transforms standard normal quantiles \eqn{z_q} to Tukey \eqn{g}-&-\eqn{h} quantiles.
#' 
#' Function [gh2z()] (compute intensive!!) transforms Tukey \eqn{g}-&-\eqn{h} quantiles to standard normal quantiles \eqn{z_q}.
#'  
#' @param z \link[base]{double} scalar or \link[base]{vector}, standard normal quantiles \eqn{z_q},
#' missingness is ***not*** allowed
#' 
#' @param q \link[base]{double} \link[base]{vector}, Tukey \eqn{g}-&-\eqn{h} quantiles \eqn{q},
#' with location parameter \eqn{A=0} and scale parameter \eqn{B=1},
#' missingness is ***not*** allowed
#' 
#' @param g,h \link[base]{double} scalars
#' 
#' @param ... parameters of function \link[rstpm2]{vuniroot}, other than `interval`, `lower` and `upper`
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
#' @importFrom rstpm2 vuniroot
#' @export
gh2z <- function(q, g = 0, h = 0, ...) {
  
  g0 <- (g == 0)
  h0 <- (h == 0)
  
  if (g0 && h0) return(q)
  
  out <- q
  
  int <- c(-8.3, 8.3)
  # identical(pnorm(8.3), 1) |> stopifnot() # even 8.29 does not work!
  
  if (!g0 && h0) { 
    # has boundary!! 
    # also has explicit form :)
    egz <- q*g + 1
    if (any(id <- (egz <= 0))) {
      out[id] <- if (g > 0) int[1L] else int[2L]
    }
    out[!id] <- log(egz[!id]) / g
    return(out)
  }
  
  if (!g0 && !h0) {
    return(vuniroot(f = \(z) expm1(g*z) * exp(h*z^2/2) - q*g, lower = int[1L], upper = int[2L], ...)[[1L]])
  }
  
  
  if (g0 && !h0) {
    return(vuniroot(f = \(z) z * exp(h*z^2/2) - q, lower = int[1L], upper = int[2L], ...)[[1L]])
  }
  
}





#' @title Derivative of [z2gh()] against `z`, on the log-scale
#' 
#' @param z \link[base]{double} \link[base]{vector}, missingness allowed
#' 
#' @param g,h \link[base]{double} scalars
#' 
#' @note
#' Nomenclature follows the parameters of ?scales::new_transform
#' 
#' @export
d_z2GH <- function(z, g, h) {
  
  hz2 <- h * z^2
  
  if (g == 0) {
    trm2 <- 1 + hz2
  } else {
    e_gz <- exp(g*z)
    trm2 <- e_gz + h * z * (e_gz - 1)/g
  }
  
  return(hz2/2 + log(trm2))
  
}










