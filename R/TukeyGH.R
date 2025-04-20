


#' @title Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description 
#' 
#' Density, distribution, quantile and simulation 
#' for Tukey \eqn{g}-&-\eqn{h} distribution with 
#' location parameter \eqn{A},
#' scale parameter \eqn{B},
#' skewness \eqn{g} and 
#' elongation \eqn{h}.
#' 
#' @param x \link[base]{double} \link[base]{vector}, quantiles \eqn{x},
#' missingness is allowed for function `fitdistrplus::fitdist()`
#' 
#' @param q \link[base]{double} \link[base]{vector}, quantiles \eqn{q}, 
#' missingness is ***not*** allowed
#' 
#' @param p \link[base]{double} \link[base]{vector}, probabilities \eqn{p},
#' missingness is ***not*** allowed
#' 
#' @param n \link[base]{integer} scalar, number of observations
#' 
#' @param log,log.p \link[base]{logical} scalar, if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' 
#' @param lower.tail \link[base]{logical} scalar, if `TRUE` (default), probabilities are \eqn{Pr(X\le x)} otherwise, \eqn{Pr(X>x)}.
#' 
#' @param A \link[base]{double} scalar, location parameter \eqn{A}, default \eqn{A=0},
#' 
#' @param B \link[base]{double} scalar, scale parameter \eqn{B>0}, default \eqn{B=1}
#' 
#' @param g \link[base]{double} scalar, skewness parameter \eqn{g}, default \eqn{g=0} (i.e., no skewness)
#' 
#' @param h \link[base]{double} scalar, elongation parameter \eqn{h\geq 0}, default \eqn{h=0} (i.e., no elongation)
#' 
#' @param ... additional parameters of function [gh2z()]
#' 
#' @returns 
#' 
#' Density function [dGH()] returns \link[base]{double} \link[base]{vector}.
#' 
#' Distribution function [pGH()] returns \link[base]{double} \link[base]{vector}.
#' 
#' Quantile function [qGH()] returns \link[base]{double} \link[base]{vector}.
#' 
#' Random generator function [rGH()] returns \link[base]{double} \link[base]{vector}.
#' 
#' @keywords internal
#' @name TukeyGH
#' @export
dGH <- function(x, A = 0, B = 1, g = 0, h = 0, log = FALSE, ...) {

  if (!(nx <- length(x))) return(double()) # ?fitdistrplus::fitdist will test len-0 `x`
  
  xok <- is.finite(x) # ?fitdistrplus::fitdist will test exceptions of x = c(0, 1, Inf, NaN, -1)
  
  z <- x
  
  if ((h < 0) || (B < 0)) { # exception handling for ?fitdistrplus::fitdist
    z[] <- NaN
    return(z)
  }
    
  z[xok] <- ((x[xok] - A)/B) |> 
    gh2z(g = g, h = h, ...)

  l <- -z^2/2 - log(2*pi)/2 - log(B) - d_z2GH(z, g = g, h = h)
  if (log) return(l) 
  return(exp(l))
  
}





# not compute-intensive
#' @rdname TukeyGH
#' @importFrom stats rnorm
#' @export
rGH <- function(n, A = 0, B = 1, g = 0, h = 0) {
  q <- n |>
    rnorm() |> 
    z2gh(g = g, h = h)
  return(A + q*B)
}


#' @rdname TukeyGH
#' @importFrom stats qnorm
#' @export
qGH <- function(p, A = 0, B = 1, g = 0, h = 0, lower.tail = TRUE, log.p = FALSE) {
  q <- p |> 
    qnorm(lower.tail = lower.tail, log.p = log.p) |>
    z2gh(g = g, h = h)
  return(A + q*B)
}


# not compute-intensive :)
#' @rdname TukeyGH
#' @importFrom stats pnorm
#' @export
pGH <- function(q, A = 0, B = 1, g = 0, h = 0, lower.tail = TRUE, log.p = FALSE, ...) {
  ((q-A)/B) |>
    gh2z(g = g, h = h, ...) |>
    pnorm(lower.tail = lower.tail, log.p = log.p)
}







