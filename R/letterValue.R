

#' @title Letter-Value Estimation of Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description
#' 
#' Letter-value based estimation (Hoaglin, 1985) of 
#' Tukey \eqn{g}-, \eqn{h}- and \eqn{g}-&-\eqn{h} distribution. 
#' All equation numbers mentioned below refer to Hoaglin (1985).
#' 
#' @param x \link[base]{double} \link[base]{vector}, one-dimensional observations
#' 
#' @param g_ \link[base]{double} \link[base]{vector}, probabilities used for estimating \eqn{g} parameter.
#' Or, use `g_ = FALSE` to implement the constraint \eqn{g=0} 
#' (i.e., an \eqn{h}-distribution is estimated). 
#' 
#' @param h_ \link[base]{double} \link[base]{vector}, probabilities used for estimating \eqn{h} parameter.
#' Or, use `h_ = FALSE` to implement the constraint \eqn{h=0}
#' (i.e., a \eqn{g}-distribution is estimated). 
#' 
#' @param halfSpread \link[base]{character} scalar, 
#' either to use `'both'` for half-spreads (default),
#' `'lower'` for half-spread, or `'upper'` for half-spread.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' Unexported function `letterV_g()` estimates parameter \eqn{g} using equation (10) for \eqn{g}-distribution
#' and the equivalent equation (31) for \eqn{g}-&-\eqn{h} distribution.
#' 
#' Unexported function `letterV_B()` estimates parameter \eqn{B} for Tukey \eqn{g}-distribution
#' (i.e., \eqn{g\neq 0}, \eqn{h=0}), using equation (8a) and (8b).
#' 
#' Unexported function `letterV_Bh_g()` estimates parameters \eqn{B} and \eqn{h} when \eqn{g\neq 0}, using equation (33).
#' 
#' Unexported function `letterV_Bh()` estimates parameters \eqn{B} and \eqn{h} for Tukey \eqn{h}-distribution,
#' i.e., when \eqn{g=0} and \eqn{h\neq 0}, using equation (26a), (26b) and (27).
#' 
#' Function [letterValue()] plays a similar role as `fitdistrplus:::start.arg.default`,
#' thus extends `fitdistrplus::fitdist` for estimating Tukey \eqn{g}-&-\eqn{h} distributions.
#' 
#' @note
#' Parameter `g_` and `h_` does not have to be truly unique; i.e., \link[base]{all.equal} elements are allowed.
#' 
#' @returns 
#' 
#' Function [letterValue()] returns a `'letterValue'` object, 
#' which is \link[base]{double} \link[base]{vector} of estimates \eqn{(\hat{A}, \hat{B}, \hat{g}, \hat{h})}
#' for a Tukey \eqn{g}-&-\eqn{h} distribution.
#' 
#' 
#' @references 
#' Hoaglin, D.C. (1985). Summarizing Shape Numerically: The \eqn{g}-and-\eqn{h} Distributions. 
#' \doi{10.1002/9781118150702.ch11}
#' 
#' @examples 
#' set.seed(77652); x = rGH(n = 1e3L, g = -.3, h = .1)
#' letterValue(x, g_ = FALSE, h_ = FALSE)
#' letterValue(x, g_ = FALSE)
#' letterValue(x, h_ = FALSE)
#' (m3 = letterValue(x))
#' 
#' library(fitdistrplus)
#' fit = fitdist(x, distr = 'GH', start = as.list.default(m3))
#' plot(fit) # fitdistrplus:::plot.fitdist
#' 
#' @name letterValue
#' @importFrom stats mad median.default quantile
#' @export
letterValue <- function(
  x,
  g_ = seq.int(from = .15, to = .25, by = .005),
  h_ = seq.int(from = .15, to = .35, by = .005),
  halfSpread = c('both', 'lower', 'upper'),
  ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  
  #### prepare parameters
  
  if (!isFALSE(g_)) {
    if (!is.double(g_) || anyNA(g_)) stop('g_ (for estimating g) must be double without missing')
    g_ <- sort.int(unique.default(g_))
    if (!length(g_)) stop('`g_` cannot be len-0')
    if (any(g_ <= 0, g_ >= .5)) stop('g_ (for estimating g) must be between 0 and .5 (not including)')
    L <- A - quantile(x, probs = g_) # lower-half-spread (LHS), the paragraph under eq.10 on p469
    U <- quantile(x, probs = 1 - g_) - A # upper-half-spread (UHS)
    ok <- (L != 0) & (U != 0) # may ==0 due to small sample size
    if (!all(ok)) g_ <- g_[ok] # else do nothing
    if (!length(g_)) stop('`g_` cannot be len-0')
  }

  if (!isFALSE(h_)) {
    if (!is.double(h_) || anyNA(h_)) stop('h_ (for estimating h) must be double without missing')
    h_ <- sort.int(unique.default(h_))
    if (!length(h_)) stop('`h_` cannot be len-0')
    if (any(h_ <= 0, h_ >= .5)) stop('h_ (for estimating h) must be between 0 and .5 (not including)')
    L <- A - quantile(x, probs = h_) # lower-half-spread (LHS), the paragraph under eq.10 on p469
    U <- quantile(x, probs = 1 - h_) - A # upper-half-spread (UHS)
    ok <- (L != 0) & (U != 0) # may ==0 due to small sample size
    if (!all(ok)) h_ <- h_[ok] # else do nothing
    if (!length(h_)) stop('`h_` cannot be len-0')
  }
  
  #### end of prepare parameters
  
  halfSpread <- match.arg(halfSpread)
  
  g <- if (isFALSE(g_)) 0 else letterV_g(x, g_ = g_)
  
  if (isFALSE(h_)) {
    if (g == 0) return(c(A = A, B = mad(x, center = A), g = 0, h = 0))
    B <- letterV_B(x, g = g, g_ = g_)
    ret <- c(A = A, B = unname(B[halfSpread]), g = g, h = 0)
    attr(ret, which = 'B') <- B
  } else {
    Bh <- if (g != 0) {
      letterV_Bh_g(x, g = g, h_ = h_) 
    } else letterV_Bh(x, h_ = h_)
    ret <- c(A = A, B = unname(Bh$B[halfSpread]), g = g, h = unname(Bh$h[halfSpread]))
    attr(ret, which = 'Bh') <- Bh
  }
  
  attr(ret, which = 'g') <- attr(g, which = 'g', exact = TRUE)
  class(ret) <- 'letterValue'
  return(ret)
  
}

#' @importFrom stats median.default .lm.fit qnorm quantile
letterV_Bh_g <- function(x, A = median.default(x), g, h_) {
  L <- A - quantile(x, probs = h_) # lower-half-spread (LHS), p469
  U <- quantile(x, probs = 1 - h_) - A # upper-half-spread (UHS)
  z <- qnorm(h_)
  regx <- z^2 / 2
  
  regyU <- log(g * U / expm1(-g*z)) # p487-eq(33)
  regyL <- log(g * L / (-expm1(g*z))) # see the equation on bottom of p486
  
  # p487, paragraph under eq(33)
  regy1 <- (regyU + regyL) / 2 # 'averaging the logrithmic results'
  regy2 <- log(g * (U+L) / (exp(-g*z) - exp(g*z))) # 'dividing the full spread with appropriate denominator', 
  # `regy2` is obtained by summing up of the numeritor and denominator of the equation on bottom of p486
  # (expm1(x) - expm1(y)) == (exp(x) - exp(y))
  # I prefer `regy2`
  regy <- regy2
  
  X <- cbind(1, regx)
  #cf1 <- unname(.lm.fit(x = X, y = cbind(regy1))$coefficients)
  #cf2 <- unname(.lm.fit(x = X, y = cbind(regy2))$coefficients)
  cf <- unname(.lm.fit(x = X, y = cbind(regy))$coefficients)
  cfL <- unname(.lm.fit(x = X, y = cbind(regyL))$coefficients)
  cfU <- unname(.lm.fit(x = X, y = cbind(regyU))$coefficients)
  
  B <- exp(c(cf[1L], cfL[1L], cfU[1L]))
  h <- pmax(0, c(cf[2L], cfL[2L], cfU[2L]))
  # slope `h` should be positive. When a negative slope is fitted, use h = 0
  names(B) <- names(h) <- c('both', 'lower', 'upper')
  
  ret <- list(B = B, h = h)
  attr(ret, which = 'Bh_g') <- data.frame(
    x = rep(regx, times = 3L), 
    y = c(regy, regyL, regyU),
    id = rep(as.character.default(1:3), each = length(h_))
    # '1' = both; '2' = lower; '3' = upper
  )
  return(ret)
}




#' @importFrom stats median.default .lm.fit qnorm quantile
letterV_Bh <- function(x, A = median.default(x), h_) {
  L <- A - quantile(x, probs = h_) # lower-half-spread (LHS), p469
  U <- quantile(x, probs = 1 - h_) - A # upper-half-spread (UHS)
  z <- qnorm(h_)
  regx <- z^2 / 2
  
  regy <- log((U+L) / (-2*z)) # p483-eq(28); all.equal(U+L, q[-id] - q[id])
  regyL <- log(-L/z) # easy to derive from p483-eq(26a)
  regyU <- log(-U/z) # easy to derive from p483-eq(26b)
  X <- cbind(1, regx)
  cf <- unname(.lm.fit(x = X, y = cbind(regy))$coefficients)
  cfL <- unname(.lm.fit(x = X, y = cbind(regyL))$coefficients)
  cfU <- unname(.lm.fit(x = X, y = cbind(regyU))$coefficients)
  
  B <- exp(c(cf[1L], cfL[1L], cfU[1L]))
  h <- pmax(0, c(cf[2L], cfL[2L], cfU[2L]))
  # slope `h` should be positive. When a negative slope is fitted, use h = 0
  names(B) <- names(h) <- c('both', 'lower', 'upper')
  ret <- list(B = B, h = h)
  attr(ret, which = 'Bh') <- data.frame(
    x = rep(regx, times = 3L), 
    y = c(regy, regyL, regyU),
    id = rep(as.character.default(1:3), each = length(h_))
    # '1' = both; '2' = lower; '3' = upper
  )
  return(ret)
}

#' @importFrom stats median.default .lm.fit qnorm quantile
letterV_B <- function(x, A = median.default(x), g, g_) {
  # when (h == 0) && (g != 0), calculate B (p469, Example, p471-eq(11a))
  q <- quantile(x, probs = c(g_, 1-g_))
  z <- qnorm(g_) # length n
  id <- seq_along(g_)
  regx <- expm1(c(g*z, -g*z))/g # p469, eq(8a-8b)
  B <- unname(.lm.fit(x = cbind(regx), y = cbind(q - A))$coefficients)
  if (B < 0) stop('estimated B < 0 ???')
  BL <- unname(.lm.fit(x = cbind(regx[id]), y = cbind(q[id] - A))$coefficients)
  if (BL < 0) stop('estimated B (lower spread) < 0 ???')
  BU <- unname(.lm.fit(x = cbind(regx[-id]), y = cbind(q[-id] - A))$coefficients)
  if (BU < 0) stop('estimated B (upper spread) < 0 ???')
  
  ret <- c(both = B, lower = BL, upper = BU)
  attr(ret, which = 'B') <- data.frame(
    x = regx, 
    y = q - A, 
    id = rep(as.character.default(1:2), each = length(g_))
    # '1' = lower; '2' = 'upper'
  )
  return(ret)
}

#' @importFrom stats median.default qnorm quantile
letterV_g <- function(x, A = median.default(x), g_) {
  # p469, equation (10) Estimating g; must use both half-spread
  L <- A - quantile(x, probs = g_) # lower-half-spread (LHS), the paragraph under eq.10 on p469
  U <- quantile(x, probs = 1 - g_) - A # upper-half-spread (UHS)
  gs <- (log(L / U) / qnorm(g_)) # p469, equation (10)
  g <- median.default(gs) 
  # p474, 2nd paragraph, the values of $g_p$ do not show regular dependence on $g_$,
  # thus we take their median as a constant-$g$ description of the skewness.
  
  attr(g, which = 'g') <- data.frame(gs = gs, p = g_) # only for plot
  return(g)
}








