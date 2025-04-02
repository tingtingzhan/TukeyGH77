

# Package OpVaR \url{https://cran.r-project.org/package=OpVaR} has been removed from CRAN.



#' @title Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description 
#' 
#' Density, distribution function, quantile function and simulation 
#' for Tukey \eqn{g}-&-\eqn{h} distribution with 
#' location parameter \eqn{A},
#' scale parameter \eqn{B},
#' skewness \eqn{g} and 
#' elongation \eqn{h}.
#' 
#' @param x,q \link[base]{double} \link[base]{vector}, quantiles
#' 
#' @param p \link[base]{double} \link[base]{vector}, probabilities
#' 
#' @param n \link[base]{integer} scalar, number of observations
#' 
#' @param log,log.p \link[base]{logical} scalar, if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' 
#' @param lower.tail \link[base]{logical} scalar, if `TRUE` (default), probabilities are \eqn{Pr(X\le x)} otherwise, \eqn{Pr(X>x)}.
#' 
#' @param A \link[base]{double} scalar, location parameter \eqn{A=0} by default
#' 
#' @param B \link[base]{double} scalar, scale parameter \eqn{B>0}. Default \eqn{B=1}
#' 
#' @param g \link[base]{double} scalar, skewness parameter \eqn{g=0} by default (i.e., no skewness)
#' 
#' @param h \link[base]{double} scalar, elongation parameter \eqn{h\geq 0}. Default \eqn{h=0} (i.e., no elongation)
#' 
# @param interval interval of standard normal quantiles, when solving from Tukey \eqn{g}-&-\eqn{h} quantiles using the vuniroot algorithm 
#' 
#' @param ... other parameters of function [vuniroot2()]
#' 
# @details
# Argument `A`, `B`, `g` and `h` will be recycled to the maximum length of the four.
#' 
#' @returns 
#' 
#' Function [dGH()] returns the density and accommodates \link[base]{vector} arguments `A`, `B`, `g` and `h`.
#' The quantiles `x` can be either \link[base]{vector} or \link[base]{matrix}.
#' This function takes about 1/5 time of `gk::dgh`.
#' 
#' Function [pGH()] returns the distribution function, only taking scalar arguments and \link[base]{vector} quantiles \eqn{q}.
#' This function takes about 1/10 time of function `gk::pgh`.
#' 
#' Function [qGH()] returns the quantile function, only taking scalar arguments and \link[base]{vector} probabilities \eqn{p}.
#' 
#' Function [rGH()] generates random deviates, only taking scalar arguments.
#' 
#' @keywords internal
#' @name TukeyGH
#' @export
dGH <- function(x, A = 0, B = 1, g = 0, h = 0, log = FALSE, ...) {
  pars <- cbind(A, B, g, h) # recycle
  return(.dGH(x = x, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], log = log, ...))
}


# not compute-intensive
#' @rdname TukeyGH
#' @importFrom stats rnorm
#' @export
rGH <- function(n, A = 0, B = 1, g = 0, h = 0) {
  n |>
    rnorm() |> 
    z2GH(A = A, B = B, g = g, h = h)
}


#' @rdname TukeyGH
#' @importFrom stats qnorm
#' @export
qGH <- function(p, A = 0, B = 1, g = 0, h = 0, lower.tail = TRUE, log.p = FALSE) {
  # only works with vector `p` and scalar `A`,`B`,`g`,`h`, for now
  z <- qnorm(p = p, lower.tail = lower.tail, log.p = log.p)
  .z2GH(z = z, A = A, B = B, g = g, h = h)
}


# not compute-intensive :)
#' @rdname TukeyGH
#' @importFrom stats pnorm
#' @export
pGH <- function(q, A = 0, B = 1, g = 0, h = 0, lower.tail = TRUE, log.p = FALSE, ...) {
  # only works with vector `q` and scalar `A`,`B`,`g`,`h`
  z <- GH2z(q = q, A = A, B = B, g = g, h = h, ...)
  pnorm(q = z, mean = 0, sd = 1, lower.tail = lower.tail, log.p = log.p)
}







