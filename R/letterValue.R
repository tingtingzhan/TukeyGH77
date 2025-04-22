

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
#' @param probs_g \link[base]{double} \link[base]{vector}, probabilities used for estimating \eqn{g} parameter.
#' Or, use `probs_g = FALSE` to implement the constraint \eqn{g=0} 
#' (i.e., an \eqn{h}-distribution is estimated). 
#' 
#' @param probs_h \link[base]{double} \link[base]{vector}, probabilities used for estimating \eqn{h} parameter.
#' Or, use `probs_h = FALSE` to implement the constraint \eqn{h=0}
#' (i.e., a \eqn{g}-distribution is estimated). 
#' 
#' @param g_select ..
#' 
#' @param true (optional) \link[base]{double} \link[base]{vector} of \eqn{(B,g,h)}, 
#' for function [letterValue_gh()] with option `g_select = 'demo'`
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
#' Parameter `probs_g` and `probs_h` does not have to be truly unique; i.e., \link[base]{all.equal} elements are allowed.
#' 
#' @returns 
#' 
#' Functions [letterValue_g()], [letterValue_h()] and [letterValue_gh()] 
#' returns a \link[base]{length}-`4L` \link[base]{double} \link[base]{vector} of
#' estimates \eqn{(\hat{A}, \hat{B}, \hat{g}, \hat{h})}.
#' 
#' @references 
#' Hoaglin, D.C. (1985). Summarizing Shape Numerically: The \eqn{g}-and-\eqn{h} Distributions. 
#' \doi{10.1002/9781118150702.ch11}
#' 
#' @keywords internal
#' @name letterValue
#' @import patchwork
#' @importFrom stats median.default
#' @export
letterValue_g <- function(
    x,
    probs_g = seq.int(from = .01, to = .3, by = .005),
    g_select = 'median', # 'median' is good enough
    ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  x <- x - A # new median = 0
  
  probs_g <- probs_g |> check_letterVal_(x = x)
  
  g <- lv_g_(x, probs_g = probs_g, g_select = g_select) 
  
  B <- lv_B_(x, g = g, probs_g = probs_g)
  
  ret <- c(A = A, B = B, g = g, h = 0)
  attr(ret, which = 'plot') <- 
    attr(g, which = 'plot', exact = TRUE) +
    attr(B, which = 'plot', exact = TRUE)
  return(ret)

}



#' @rdname letterValue
#' @import patchwork
#' @importFrom stats median.default
#' @export
letterValue_h <- function(
    x,
    probs_h = seq.int(from = .1, to = .3, by = .005),
    ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  x <- x - A # new median = 0
  
  probs_h <- probs_h |> check_letterVal_(x = x)
  
  Bh <- lv_Bh_(x, g = 0, probs_h = probs_h)
  
  ret <- c(A = A, B = unname(Bh['B']), g = 0, h = unname(Bh['h']))
  attr(ret, which = 'plot') <- attr(Bh, which = 'plot', exact = TRUE)
  return(ret)
  
}


#' @rdname letterValue
#' @import patchwork
#' @importFrom stats median.default
#' @export
letterValue_gh <- function(
  x,
  probs_g = seq.int(from = .01, to = .3, by = .005),
  probs_h = seq.int(from = .1, to = .3, by = .005),
  g_select = 'h_optim',
  #true, 
  ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  x <- x - A # new median = 0
  
  probs_g <- probs_g |> check_letterVal_(x = x)
  probs_h <- probs_h |> check_letterVal_(x = x)
  
  g <- lv_g_(x, probs_g = probs_g, probs_h = probs_h, g_select = g_select, ...)
  
  if (g_select == 'demo') {
    
    #if (missing(true)) stop('must provide true `B`, `g` and `h`')
    
    demo_plot <- mapply(FUN = \(g, nm, pane) {
      names(g) <- nm
      (lv_Bh_(x, g = g, probs_h = probs_h, ...) |> 
          attr(which = 'plot', exact = TRUE)) +
        labs(title = sprintf(fmt = '(%s). %s', pane, nm))
    }, g = g, nm = names(g), pane = LETTERS[seq_along(g)], SIMPLIFY = FALSE) |>
      Reduce(f = `+`)
    
    demo_plot <- demo_plot + 
      (attr(g, which = 'plot', exact = 'TRUE') + 
         labs(title = sprintf(fmt = '(%s). Estimated g', LETTERS[length(g) + 1L])))
    
    return(demo_plot + plot_layout(ncol = 3L))
    
  }
  
  Bh <- lv_Bh_(x, g = g, probs_h = probs_h)
  
  ret <- c(A = A, B = unname(Bh['B']), g = unname(g), h = unname(Bh['h']))
  attr(ret, which = 'plot') <- attr(g, which = 'plot', exact = TRUE) +
    attr(Bh, which = 'plot', exact = TRUE)
  return(ret)
  
}



#' @importFrom stats quantile
check_letterVal_ <- function(x, probs) {
  if (!is.double(probs) || anyNA(probs)) stop('`probs` must be double without missing')
  probs <- probs |> unique.default() |> sort.int()
  if (!length(probs)) stop('`probs` cannot be len-0')
  if (any(probs <= 0, probs >= .5)) stop('`probs` must be between 0 and .5 (not including)')
  lhs <- - quantile(x, probs = probs) # lower-half-spread (LHS), the paragraph under eq.10 on p469
  uhs <- quantile(x, probs = 1 - probs) # upper-half-spread (UHS)
  gok <- (lhs != 0) & (uhs != 0) # may ==0 due to small sample size
  if (!any(gok)) stop('data degeneration')
  probs <- probs[gok]
  attr(probs, which = 'lhs') <- lhs[gok]
  attr(probs, which = 'uhs') <- uhs[gok]
  return(probs)
}









#' @importFrom geomtextpath geom_textabline geom_textsmooth
#' @importFrom ggplot2 ggplot geom_point geom_smooth geom_abline scale_color_discrete scale_fill_discrete scale_x_continuous sec_axis xlim labs theme
#' @importFrom grid unit
#' @importFrom stats .lm.fit qnorm quantile
#' @importFrom scales label_percent
#' @importFrom ggplot2 xlim
lv_Bh_ <- function(
    x, 
    g, 
    probs_h,
    do.plot = TRUE,
    true,
    ...
) {
  
  lhs <- attr(probs_h, which = 'lhs', exact = TRUE)
  uhs <- attr(probs_h, which = 'uhs', exact = TRUE)
  z <- qnorm(probs_h)
  regx <- z^2 / 2
  
  if (g == 0) {
    
    regy <- log((uhs+lhs) / (-2*z)) # p483-eq(28); all.equal(uhs+lhs, q[-id] - q[id])
    regyL <- log(-lhs/z) # easy to derive from p483-eq(26a)
    regyU <- log(-uhs/z) # easy to derive from p483-eq(26b)
  
  } else {
    
    regyU <- log(g * uhs / expm1(-g*z)) # p487-eq(33)
    regyL <- log(g * lhs / (-expm1(g*z))) # see the equation on bottom of p486
    
    # p487, paragraph under eq(33)
    regy1 <- (regyU + regyL) / 2 # 'averaging the logrithmic results'
    regy2 <- log(g * (uhs+lhs) / (exp(-g*z) - exp(g*z))) # 'dividing the full spread with appropriate denominator', 
    # `regy2` is obtained by summing up of the numeritor and denominator of the equation on bottom of p486
    # (expm1(x) - expm1(y)) == (exp(x) - exp(y))
    # I prefer `regy2`
    regy <- regy2
    
  }
  
  X <- cbind(1, regx)
  cf <- unname(.lm.fit(x = X, y = cbind(regy))$coefficients)
  cfL <- unname(.lm.fit(x = X, y = cbind(regyL))$coefficients)
  cfU <- unname(.lm.fit(x = X, y = cbind(regyU))$coefficients)
  
  # c('both', 'lower', 'upper')
  B <- exp(c(cf[1L], cfL[1L], cfU[1L]))
  h <- c(cf[2L], cfL[2L], cfU[2L])
  # slope `h` should be positive. When a negative slope is fitted, use h = 0

  ret <- c(B = B[1L], h = pmax(0, h[1L]))
  attr(ret, which = 'errB') <- (cfL[1L] - cfU[1L])^2
  attr(ret, which = 'errh') <- (cfL[2L] - cfU[2L])^2
  if (!do.plot) return(ret) # when doing `g`-selection
  
  hs <- 1:3 |> # c('both', 'lower', 'upper')
    as.character() |>
    rep(each = length(probs_h))
  est <- sprintf(fmt = '%s\nB=%.3f\nh=%.3f', c('Both', 'Lower', 'Upper'), B, h)
  mp <- aes(
    x = rep(regx, times = 3L), 
    y = c(regy, regyL, regyU), 
    #label = est |> rep(each = length(probs_h)),
    color = hs, fill = hs)
  
  mp_smooth <- aes(
    x = rep(regx, times = 3L), 
    y = c(regy, regyL, regyU), 
    label = sprintf(fmt = '%s: B=%.3f, h=%.3f', c('Both', 'Lower', 'Upper'), B, h) |> rep(each = length(probs_h)),
    color = hs)
  
  attr(ret, which = 'plot') <- ggplot() + 
    geom_point(mapping = mp, alpha = .2, show.legend = FALSE) + 
    #geom_smooth(mapping = mp, formula = y ~ x, alpha = .2, linewidth = .5, linetype = 2L, method = 'lm', se = FALSE) +
    geom_textsmooth(
      mapping = mp_smooth, formula = y ~ x, linewidth = .5, linetype = 2L, method = 'lm', show.legend = FALSE
    ) +
    # (if (!missing(true)) geom_abline(intercept = log(true['B']), slope = true['h'], colour = 'grey40', linetype = 2L)) +
    (if (!missing(true)) geom_textabline(
      intercept = log(true['B']), slope = true['h'], 
      label = sprintf(fmt = 'true: B = %.2g, h = %.2g', true['B'], true['h']), 
      colour = 'grey40', linetype = 2L)) + 
    #scale_color_discrete(name = NULL, breaks = 1:3, labels = est) +
    #scale_fill_discrete(name = NULL, breaks = 1:3, labels = est) +
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
        'Constrained at g = 0'
      } else sprintf(fmt = '%s g = %.3f', names(g), g)
    ) + 
    theme(
      legend.position = 'inside',
      legend.position.inside = c(.8, .3),
      legend.key.spacing.y = unit(.01, units = 'npc')
    )
  
  return(ret)
}













#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_x_continuous labs annotate
#' @importFrom latex2exp TeX
#' @importFrom scales pal_hue
#' @importFrom stats .lm.fit qnorm quantile
lv_B_ <- function(x, g, probs_g) {
  # when (h == 0) && (g != 0), calculate B (p469, Example, p471-eq(11a))
  q <- quantile(x, probs = c(probs_g, 1-probs_g)) - 0
  z <- qnorm(probs_g) # length n
  id <- seq_along(probs_g)
  regx <- expm1(c(g*z, -g*z))/g # p469, eq(8a-8b)
  B <- unname(.lm.fit(x = cbind(regx), y = cbind(q))$coefficients)
  if (B < 0) stop('estimated B < 0')
  BL <- unname(.lm.fit(x = cbind(regx[id]), y = cbind(q[id]))$coefficients)
  if (BL < 0) stop('estimated B (lower spread) < 0')
  BU <- unname(.lm.fit(x = cbind(regx[-id]), y = cbind(q[-id]))$coefficients)
  if (BU < 0) stop('estimated B (upper spread) < 0')
  
  hs <- c('lower', 'upper') |>
    rep(each = length(probs_g))
  mp <- aes(x = regx, y = q, color = hs, fill = hs)
    
  attr(B, which = 'plot') <- ggplot() + 
    geom_point(mapping = mp, alpha = .2, show.legend = FALSE) + 
    geom_smooth(mapping = mp, formula = y ~ x, alpha = .2, linewidth = .5, linetype = 2L, method = 'lm', se = FALSE, show.legend = FALSE) +
    geom_smooth(mapping = aes(x = regx, y = q), formula = y ~ x, alpha = .2, colour = 'grey90', linewidth = .3, linetype = 3L, method = 'lm', se = FALSE, show.legend = FALSE) +
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
      #label = sprintf(fmt = 'Slope$\\hat{B} = %.3f$', c(BL, BU, B)) |> TeX() # why warning??
      label = sprintf(fmt = 'Slope\nB=%.3f', c(BL, BU, B))
    ) +
    labs(
      y = '$t_p$' |> TeX(), 
      caption = sprintf(fmt = '%s g = %.3f; h = 0', names(g), g)
    )
  
  return(B)
}







#' @importFrom stats optimize qnorm quantile
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_y_continuous dup_axis guide_axis labs
#' @importFrom latex2exp TeX
#' @importFrom scales label_percent
lv_g_ <- function(
    x, 
    probs_g, 
    probs_h,
    g_select = c('median', 'B_optim', 'h_optim', 'demo'),
    true,
    ...
) {
  
  lhs <- attr(probs_g, which = 'lhs', exact = TRUE)
  uhs <- attr(probs_g, which = 'uhs', exact = TRUE)
  # p469, equation (10) Estimating g; must use both half-spread
  gs <- (log(lhs / uhs) / qnorm(probs_g)) # p469, equation (10)
  
  g_q4 <- quantile(gs, probs = c(.25, .5, .75))
  
  g_select <- match.arg(g_select)
  
  if (g_select %in% c('median', 'demo')) {
    g_median_ <- unname(g_q4[2L])
    # p474, 2nd paragraph, the values of $g_p$ do not show regular dependence on $probs_g$,
    # thus we take their median as a constant-$g$ description of the skewness.
  }
  
  if (g_select %in% c('B_optim', 'demo')) {
    g_B_ <- optimize(f = \(g) {
      lv_Bh_(x = x, g = g, probs_h = probs_h, do.plot = FALSE) |>
        attr(which = 'errB', exact = TRUE)
    }, lower = g_q4[1L], upper = g_q4[3L])$minimum
  }
  
  if (g_select %in% c('h_optim', 'demo')) {
    g_h_ <- optimize(f = \(g) {
      lv_Bh_(x = x, g = g, probs_h = probs_h, do.plot = FALSE) |>
        attr(which = 'errh', exact = TRUE)
    }, lower = g_q4[1L], upper = g_q4[3L])$minimum
  }
  
  g <- switch(g_select, median = {
    c(median = g_median_)
  }, B_optim = {
    c(B_optim = g_B_)
  }, h_optim = {
    c(h_optim = g_h_)
  }, demo = {
    c(median = g_median_, B_optim = g_B_, h_optim = g_h_, true = unname(true['g']))
  })
  
  attr(g, which = 'plot') <- ggplot() + 
    geom_point(mapping = aes(x = probs_g, y = gs)) + 
    scale_x_continuous(name = 'p', labels = label_percent()) +
    scale_y_continuous(
      name = '$\\hat{g}_p$' |> TeX(),
      sec.axis = dup_axis(
        name = NULL, 
        breaks = g, 
        labels = sprintf(fmt = '%s\n%.3f', names(g), g),
        guide = guide_axis(angle = switch(g_select, demo = 45, 90))
      )
    )
  
  return(g)
}








