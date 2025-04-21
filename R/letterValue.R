

#' @title Letter-Value Estimation of Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description
#' 
#' Letter-value based estimation (Hoaglin, 1985) of 
#' Tukey \eqn{g}-, \eqn{h}- and \eqn{g}-&-\eqn{h} distribution. 
#' All equation numbers mentioned below refer to Hoaglin (1985).
#' 
#' @param x \link[base]{double} \link[base]{vector}, observations \eqn{x},
#' missingness ***not*** allowed.
#' 
#' @param g_ \link[base]{double} \link[base]{vector}, probabilities used for estimating \eqn{g} parameter.
#' Or, use `g_ = FALSE` to implement the constraint \eqn{g=0} 
#' (i.e., an \eqn{h}-distribution is estimated). 
#' 
#' @param h_ \link[base]{double} \link[base]{vector}, probabilities used for estimating \eqn{h} parameter.
#' Or, use `h_ = FALSE` to implement the constraint \eqn{h=0}
#' (i.e., a \eqn{g}-distribution is estimated). 
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' Unexported function `lv_g_()` estimates parameter \eqn{g} using equation (10) for \eqn{g}-distribution
#' and the equivalent equation (31) for \eqn{g}-&-\eqn{h} distribution.
#' 
#' Unexported function `lv_B_()` estimates parameter \eqn{B} for Tukey \eqn{g}-distribution
#' (i.e., \eqn{g\neq 0}, \eqn{h=0}), using equation (8a) and (8b).
#' 
#' Unexported function `lv_Bh_()` estimates parameters \eqn{B} and \eqn{h} 
#' with given \eqn{g\neq 0}, using equation (33), 
#' or with given \eqn{g=0}, using equation (26a), (26b) and (27).
#' 
#' Unexported function `letterV_Bh()` estimates parameters \eqn{B} and \eqn{h} for Tukey \eqn{h}-distribution,
#' i.e., with given \eqn{g=0}, using equation (26a), (26b) and (27).
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
#' @keywords internal
#' @name letterValue
#' @import patchwork
#' @importFrom stats mad median.default quantile
#' @export
letterValue <- function(
  x,
  g_ = seq.int(from = .01, to = .3, by = .005),
  h_ = seq.int(from = .1, to = .3, by = .005),
  ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  x <- x - A # new median = 0
  
  if (!isFALSE(g_)) {
    if (!is.double(g_) || anyNA(g_)) stop('g_ (for estimating g) must be double without missing')
    g_ <- g_ |> unique.default() |> sort.int()
    if (!length(g_)) stop('`g_` cannot be len-0')
    if (any(g_ <= 0, g_ >= .5)) stop('g_ (for estimating g) must be between 0 and .5 (not including)')
    gL <- 0 - quantile(x, probs = g_) # lower-half-spread (LHS), the paragraph under eq.10 on p469
    gU <- quantile(x, probs = 1 - g_) - 0 # upper-half-spread (UHS)
    gok <- (gL != 0) & (gU != 0) # may ==0 due to small sample size
    if (!any(gok)) stop('data degeneration')
    if (!all(gok)) {
      g_ <- g_[gok]
      gL <- gL[gok]
      gU <- gU[gok]
    } # else do nothing
  }

  if (!isFALSE(h_)) {
    if (!is.double(h_) || anyNA(h_)) stop('h_ (for estimating h) must be double without missing')
    h_ <- h_ |> unique.default() |> sort.int()
    if (!length(h_)) stop('`h_` cannot be len-0')
    if (any(h_ <= 0, h_ >= .5)) stop('h_ (for estimating h) must be between 0 and .5 (not including)')
    hL <- 0 - quantile(x, probs = h_) # lower-half-spread (LHS), the paragraph under eq.10 on p469
    hU <- quantile(x, probs = 1 - h_) - 0 # upper-half-spread (UHS)
    hok <- (hL != 0) & (hU != 0) # may ==0 due to small sample size
    if (!any(hok)) stop('data degeneration')
    if (!all(hok)) {
      h_ <- h_[hok] 
      hL <- hL[hok]
      hU <- hU[hok]
    } # else do nothing
  }
  
  g <- if (isFALSE(g_)) 0 else lv_g_(x, g_ = g_, gL = gL, gU = gU)
  
  if (isFALSE(h_)) {
    
    if (g == 0) {
      return(c(A = A, B = mad(x, center = 0), g = 0, h = 0))
    }
    
    B <- lv_B_(x, g = g, g_ = g_)
    ret <- c(A = A, B = B, g = g, h = 0)
    attr(ret, which = 'plot') <- 
      attr(g, which = 'plot', exact = TRUE) +
      attr(B, which = 'plot', exact = TRUE)
    return(ret)
    
  } 
    
  Bh <- lv_Bh_(x, g = g, h_ = h_, hL = hL, hU = hU) 
  ret <- c(A = A, B = unname(Bh['B']), g = g, h = unname(Bh['h']))
  
  attr(ret, which = 'plot') <- if (g == 0) {
    attr(Bh, which = 'plot', exact = TRUE)
  } else {
    attr(g, which = 'plot', exact = TRUE) +
    attr(Bh, which = 'plot', exact = TRUE)
  }
  
  return(ret)
  
}



#' @importFrom ggplot2 ggplot geom_point geom_smooth scale_color_discrete scale_fill_discrete scale_x_continuous sec_axis xlim labs theme
#' @importFrom grid unit
#' @importFrom stats .lm.fit qnorm quantile
#' @importFrom scales label_percent
#' @importFrom ggplot2 xlim
lv_Bh_ <- function(x, g, h_, hL, hU) {
  
  z <- qnorm(h_)
  regx <- z^2 / 2
  
  if (g == 0) {
    
    regy <- log((hU+hL) / (-2*z)) # p483-eq(28); all.equal(hU+hL, q[-id] - q[id])
    regyL <- log(-hL/z) # easy to derive from p483-eq(26a)
    regyU <- log(-hU/z) # easy to derive from p483-eq(26b)
    
    subtitle <- c( # double check the code vs. label!
      'Combined: $log(-(x_{1-p} - x_p)/(2z_p))',
      'Lower: $log(x_p/z_p)', 
      'Upper: $log(-x_{1-p})/z_p)'
    )
    
  } else {
    
    regyU <- log(g * hU / expm1(-g*z)) # p487-eq(33)
    regyL <- log(g * hL / (-expm1(g*z))) # see the equation on bottom of p486
    
    # p487, paragraph under eq(33)
    regy1 <- (regyU + regyL) / 2 # 'averaging the logrithmic results'
    regy2 <- log(g * (hU+hL) / (exp(-g*z) - exp(g*z))) # 'dividing the full spread with appropriate denominator', 
    # `regy2` is obtained by summing up of the numeritor and denominator of the equation on bottom of p486
    # (expm1(x) - expm1(y)) == (exp(x) - exp(y))
    # I prefer `regy2`
    regy <- regy2
    
    subtitle <- c( # double check the code vs. label!
      'Lower: $log(-g(-x_p)/(exp(gz_p) - 1))',
      'Upper: $log(g(x_{1-p})/(exp(-gz_p) - 1))',
      # 'Combo v1: averaging response of Lower & Upper',
      # 'Combo v2: $log(g(x_{1-p} - x_p)/(exp(-gz_p)-exp(gz)))'
      'Combo: $log(g(x_{1-p} - x_p)/(exp(-gz_p)-exp(gz)))'
    )
  }
  
  X <- cbind(1, regx)
  cf <- unname(.lm.fit(x = X, y = cbind(regy))$coefficients)
  cfL <- unname(.lm.fit(x = X, y = cbind(regyL))$coefficients)
  cfU <- unname(.lm.fit(x = X, y = cbind(regyU))$coefficients)
  
  # c('both', 'lower', 'upper')
  B <- exp(c(cf[1L], cfL[1L], cfU[1L]))
  h <- pmax(0, c(cf[2L], cfL[2L], cfU[2L]))
  # slope `h` should be positive. When a negative slope is fitted, use h = 0

  ret <- c(B = B[1L], h = h[1L])
  
  #nms <- sprintf(fmt = '%s\nB\u0302=%.3f, h\u0302=%.3f', c(
  #  ### 'Combo (average log)', 
  #  ### 'Combo (full-spread)', 
  #  'Combo',
  #  'Lower', 'Upper'), Bh$B, Bh$h)
  
  hs <- 1:3 |> # c('both', 'lower', 'upper')
    as.character() |>
    rep(each = length(h_))
  est <- sprintf(fmt = 'B=%.3f\nh=%.3f', B, h)
  mp <- aes(
    x = rep(regx, times = 3L), 
    y = c(regy, regyL, regyU), 
    color = hs, fill = hs)
  
  attr(ret, which = 'plot') <- ggplot() + 
    geom_point(mapping = mp, alpha = .2) + 
    geom_smooth(mapping = mp, formula = y ~ x, alpha = .2, linewidth = .5, linetype = 2L, method = 'lm') +
    scale_color_discrete(name = NULL, breaks = 1:3, labels = est) +
    scale_fill_discrete(name = NULL, breaks = 1:3, labels = est) +
    scale_x_continuous(
      name = '$z_p^2/2$' |> TeX(),
      sec.axis = sec_axis(
        name = 'p', 
        transform = ~ pnorm(sqrt(. * 2), lower.tail = FALSE),
        labels = label_percent()
      )
    ) +
    labs(
      y = NULL,
      caption = if (g == 0) {
        'Given g = 0\nCombined: p483, eq(27)\nUpper & Lower: p483, eq(26a-26b)'
      } else sprintf(fmt = 'Given g = %.3f\nUpper: p487,eq(33)\nLower: see p486,bottom', g)
    ) + 
    theme(
      legend.position = 'inside',
      legend.position.inside = c(.8, .2),
      legend.key.spacing.y = unit(.01, units = 'npc')
    )
  
  return(ret)
}













#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_x_continuous labs annotate
#' @importFrom latex2exp TeX
#' @importFrom scales pal_hue
#' @importFrom stats .lm.fit qnorm quantile
lv_B_ <- function(x, g, g_) {
  # when (h == 0) && (g != 0), calculate B (p469, Example, p471-eq(11a))
  q <- quantile(x, probs = c(g_, 1-g_)) - 0
  z <- qnorm(g_) # length n
  id <- seq_along(g_)
  regx <- expm1(c(g*z, -g*z))/g # p469, eq(8a-8b)
  B <- unname(.lm.fit(x = cbind(regx), y = cbind(q))$coefficients)
  if (B < 0) stop('estimated B < 0')
  BL <- unname(.lm.fit(x = cbind(regx[id]), y = cbind(q[id]))$coefficients)
  if (BL < 0) stop('estimated B (lower spread) < 0')
  BU <- unname(.lm.fit(x = cbind(regx[-id]), y = cbind(q[-id]))$coefficients)
  if (BU < 0) stop('estimated B (upper spread) < 0')
  
  hs <- c('lower', 'upper') |>
    rep(each = length(g_))
  mp <- aes(x = regx, y = q, color = hs, fill = hs)
    
  attr(B, which = 'plot') <- ggplot() + 
    geom_point(mapping = mp, alpha = .2, show.legend = FALSE) + 
    geom_smooth(mapping = mp, formula = y ~ x, alpha = .2, linewidth = .5, linetype = 2L, method = 'lm', show.legend = FALSE) +
    geom_smooth(mapping = aes(x = regx, y = q), formula = y ~ x, alpha = .2, colour = 'grey90', linewidth = .3, linetype = 3L, method = 'lm', show.legend = FALSE) +
    scale_x_continuous(
      name = '$(e^{gz_p} - 1)/g$' |> TeX(),
      sec.axis = sec_axis(
        name = 'p', 
        transform = ~ pnorm(log(.*g+1)/g),
        breaks = c(.05, .1, .25, .5, .75, .9, .95),
        labels = label_percent()
      )
    ) +
    annotate(
      geom = 'text', 
      x = c(mean(regx[id]), mean(regx[-id]), mean(regx)), y = c(mean(q[id]), mean(q[-id]), mean(q)), 
      colour = c(pal_hue()(2L), 'grey40'),
      fontface = c('plain', 'plain', 'bold'),
      #label = sprintf(fmt = '$\\hat{B} = %.3f$', c(BL, BU, B)) |> TeX() # why warning??
      label = sprintf(fmt = 'Slope\nB=%.3f', c(BL, BU, B))
    ) +
    labs(
      y = '$t_p$' |> TeX(), 
      caption = sprintf(fmt = 'Given g = %.3f; h = 0\nLower Half Spread: p469, eq(8a)\nUpper Half Spread: p469, eq(8b)', g)
    )
  
  return(B)
}







#' @importFrom stats qnorm quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_x_continuous scale_y_continuous dup_axis guide_axis labs
#' @importFrom latex2exp TeX
#' @importFrom scales label_percent
lv_g_ <- function(x, g_, gL, gU) {
  # p469, equation (10) Estimating g; must use both half-spread
  gs <- (log(gL / gU) / qnorm(g_)) # p469, equation (10)
  g <- median.default(gs) 
  # p474, 2nd paragraph, the values of $g_p$ do not show regular dependence on $g_$,
  # thus we take their median as a constant-$g$ description of the skewness.
  
  attr(g, which = 'plot') <- ggplot() + 
    geom_point(mapping = aes(x = g_, y = gs)) + 
    geom_hline(yintercept = g, linetype = 2L, alpha = .2) +
    scale_x_continuous(name = 'p', labels = label_percent()) +
    scale_y_continuous(
      name = '$\\hat{g}_p$' |> TeX(),
      sec.axis = dup_axis(
        name = NULL, 
        breaks = g, 
        labels = \(i) sprintf(fmt = 'median\n%.3f', i),
        guide = guide_axis(angle = 90)
      )
    ) +
    labs(caption = 'p469, eq(10)')
  
  return(g)
}








