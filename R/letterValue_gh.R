

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
#' @param A_select ..
#' 
#' @param g_select ..
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
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
letterValue_gh <- function(
  x,
  probs_g = seq.int(from = .01, to = .3, by = .005),
  probs_h = seq.int(from = .1, to = .3, by = .005),
  A_select = c('optim', 'median', 'demo'),
  g_select = c('h.optim', 'median', 'B.optim', 'demo'),
  #do.plot = FALSE,
  ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  
  probs_g <- probs_g |> check_letterVal_(x = x)
  probs_h <- probs_h |> check_letterVal_(x = x)
  
  A_select <- match.arg(A_select)
  if (A_select == 'demo') {
    p <- mapply(FUN = \(s, pane) {
      p <- lv_g_A(x = x, probs_g = probs_g, A_select = s) |>
        attr(which = 'plot', exact = TRUE)
      p$labels$title <- sprintf(fmt = '(%s)', pane) |> 
        paste(p$labels$title)
      return(p)
    }, s = c('median', 'optim'), pane = LETTERS[seq_len(2L)], SIMPLIFY = FALSE) |>
      Reduce(f = `+`)
    return(p)
  }
  A <- lv_g_A(x = x, probs_g = probs_g, A_select = A_select, ...)
  
  g_select <- match.arg(g_select)
  g <- lv_gh_g(x, A = A, probs_g = probs_g, probs_h = probs_h, g_select = g_select, ...)
  if (g_select == 'demo') {
    p_g <- g |> 
      attr(which = 'plot', exact = TRUE)
    p_Bh <- c('median', 'B.optim', 'h.optim') |>
      lapply(FUN = \(s) {
        g <- lv_gh_g(x, A = A, probs_g = probs_g, probs_h = probs_h, g_select = s, ...) 
        lv_Bh_(x = x, A = A, g = g, probs_g = probs_g, probs_h = probs_h) |>
          attr(which = 'plot', exact = TRUE)
      })
    p <- p_g |>
      list() |>
      c(p_Bh) |>
      add_pane() |>
      Reduce(f = `+`)
    return(p)
  }
  
  Bh <- lv_Bh_(x, A = A, g = g, probs_h = probs_h)
  
  ret <- c(A = A, B = unname(Bh['B']), g = unname(g), h = unname(Bh['h']))
  #if (do.plot) {
  #  attr(ret, which = 'plot') <- attr(g, which = 'plot', exact = TRUE) + attr(Bh, which = 'plot', exact = TRUE)
  #}
  return(ret)
  
}



#' @importFrom stats quantile
check_letterVal_ <- function(x, probs) {
  if (!is.double(probs) || anyNA(probs)) stop('`probs` must be double without missing')
  probs <- probs |> unique.default() |> sort.int()
  if (!length(probs)) stop('`probs` cannot be len-0')
  if (any(probs <= 0, probs >= .5)) stop('`probs` must be between 0 and .5 (not including)')
  
  # bound LHS from 40th-percentile
  # bound UHS from 60th-percentile
  # very conservative
  lhs <- quantile(x, probs = .4) - quantile(x, probs = probs) # lower-half-spread (LHS), the paragraph under eq.10 on p469
  uhs <- quantile(x, probs = 1 - probs) - quantile(x, probs = .6) # upper-half-spread (UHS)
  
  gok <- (lhs > 0) & (uhs > 0) # may due to small sample size
  if (!any(gok)) stop('data degeneration')
  
  probs <- probs[gok]
  return(probs)
}









#' @importFrom geomtextpath geom_textabline geom_textsmooth
#' @importFrom ggplot2 ggplot geom_point geom_abline scale_x_continuous sec_axis labs
#' @importFrom stats .lm.fit qnorm quantile
#' @importFrom scales label_percent
lv_Bh_ <- function(
    x, A, g, 
    probs_h,
    ...
) {
  
  tL <- quantile(x, probs = probs_h)
  tU <- quantile(x, probs = 1 - probs_h)
  z <- qnorm(probs_h)
  n <- length(probs_h)
  
  r_x <- z^2 / 2
  r_yL <- log(g * (tL-A) / expm1(g*z))
  r_yU <- log(g * (tU-A) / expm1(-g*z))
  r_y <- log(g * (tU-tL) / (exp(-g*z) - exp(g*z)))
    
  X <- cbind(1, r_x)
  cf <- unname(.lm.fit(x = X, y = cbind(r_y))$coefficients)
  cfL <- unname(.lm.fit(x = X, y = cbind(r_yL))$coefficients)
  cfU <- unname(.lm.fit(x = X, y = cbind(r_yU))$coefficients)
  
  # c('both', 'lower', 'upper')
  B <- exp(cf[1L])
  h <- cf[2L]
  # slope `h` should be positive. When a negative slope is fitted, use h = 0

  ret <- c(B = B, h = pmax(0, h))
  
  aes_d_ <- data.frame(
    x = r_x,
    y = c(r_yL, r_yU),
    label = sprintf(fmt = '%s: B=%.3f, h=%.3f', c('Lower', 'Upper'), exp(c(cfL[1L], cfU[1L])), c(cfL[2L], cfU[2L])) |> 
      rep(each = n)
  )
  
  attr(ret, which = 'plot') <- ggplot() + 
    
    geom_point(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, color = .data$label), alpha = .2, show.legend = FALSE) + 
    
    geom_textsmooth(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, label = .data$label, color = .data$label), formula = y ~ x, linetype = 2L, method = 'lm', show.legend = FALSE) +
    
    geom_textsmooth(mapping = aes(x = r_x, y = r_y), 
                    label = sprintf(fmt = 'B=%.3f, h=%.3f', B, h), 
                    color = 'grey40', formula = y ~ x, linetype = 2L, method = 'lm', show.legend = FALSE) +
    
    scale_x_continuous(
      name = '$z_p^2/2$' |> TeX(),
      sec.axis = sec_axis(
        name = 'p', 
        transform = ~ pnorm(sqrt(. * 2), lower.tail = FALSE),
        labels = label_percent()
      )
    ) +
    
    labs(
      subtitle = sprintf(fmt = '$\\hat{g}_{%s} = %.3f$', names(g), g) |> TeX(),
      y = NULL
    )
  
  return(ret)
}


















#' @importFrom stats optimize qnorm quantile
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_y_continuous dup_axis guide_axis labs
#' @importFrom latex2exp TeX
#' @importFrom scales label_percent
lv_gh_g <- function(
    x, A,
    probs_g, 
    probs_h,
    g_select = c('median', 'B.optim', 'h.optim', 'demo'),
    ...
) {
  
  gs <- log((A - quantile(x, probs = probs_g)) / (quantile(x, probs = 1 - probs_g) - A)) / qnorm(probs_g)
  g_intv_ <- quantile(gs, probs = c(.4, .6))
  
  g_select <- match.arg(g_select)
  
  if (g_select %in% c('median', 'demo')) {
    g_median_ <- median.default(gs)
  }
  
  tL <- quantile(x, probs = probs_h)
  tU <- quantile(x, probs = 1 - probs_h)
  z <- qnorm(probs_h)
  foo <- \(g) {
    r_x <- z^2 / 2
    r_yL <- log(g * (tL-A) / expm1(g*z))
    r_yU <- log(g * (tU-A) / expm1(-g*z))
    r_X <- cbind(1, r_x)
    cfL <- unname(.lm.fit(x = r_X, y = cbind(r_yL))$coefficients)
    cfU <- unname(.lm.fit(x = r_X, y = cbind(r_yU))$coefficients)
    ret <- (cfL - cfU)^2
    names(ret) <- c('errB', 'errh')
    return(ret)
  }

  if (g_select %in% c('B.optim', 'demo')) {
    g_B_ <- optimize(f = \(g) foo(g)['errB'], interval = g_intv_)$minimum
  }
  
  if (g_select %in% c('h.optim', 'demo')) {
    g_h_ <- optimize(f = \(g) foo(g)['errh'], interval = g_intv_)$minimum
  }
  
  g <- switch(g_select, median = {
    c(median = g_median_)
  }, B.optim = {
    c(B.optim = g_B_)
  }, h.optim = {
    c(h.optim = g_h_)
  }, demo = {
    c(median = g_median_, B.optim = g_B_, h.optim = g_h_)
  })
  
  attr(g, which = 'plot') <- attr(A, which = 'plot', exact = TRUE) +
    scale_y_continuous(
      sec.axis = dup_axis(
        name = NULL, 
        breaks = g, 
        labels = sprintf(fmt = '%s\n%.3f', names(g), g),
        guide = guide_axis(angle = switch(g_select, demo = 45, 90))
      )
    )
  
  return(g)
}






