

#' @rdname letterValue
#' @import patchwork
#' @importFrom stats median.default
#' @export
letterValue_g <- function(
    x,
    probs_g = seq.int(from = .01, to = .3, by = .005),
    A_select = c('optim', 'median', 'demo'),
    g_select = c('optim', 'median', 'demo'),
    #do.plot = FALSE,
    ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  
  probs_g <- probs_g |> check_letterVal_(x = x)
  
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
  g <- lv_g_g(x, A = A, probs_g = probs_g, g_select = g_select, ...) 
  if (g_select == 'demo') {
    p_g <- g |> 
      attr(which = 'plot', exact = TRUE)
    p_B <- c('median', 'optim') |>
      lapply(FUN = \(s) {
        g <- lv_g_g(x, A = A, probs_g = probs_g, g_select = s, ...) 
        lv_B_(x = x, A = A, g = g, probs_g = probs_g) |>
          attr(which = 'plot', exact = TRUE)
      })
    p <- p_g |>
      list() |>
      c(p_B) |>
      add_pane() |>
      Reduce(f = `+`)
    return(p)
  }
  
  B <- lv_B_(x, A = A, g = g, probs_g = probs_g)
  
  ret <- c(A = A, B = B, g = unname(g), h = 0)
  #if (do.plot) {
  #  attr(ret, which = 'plot') <- list(
  #    attr(g, which = 'plot', exact = TRUE),
  #    attr(B, which = 'plot', exact = TRUE)
  #  ) |> 
  #    add_pane() |> 
  #    Reduce(f = `+`)
  #}
  return(ret)
  
}



add_pane <- function(p, pane = LETTERS[seq_along(p)], ...) {
  mapply(FUN = \(p, pane) {
    p$labels$title <- sprintf(fmt = '(%s)', pane) |> 
      paste(p$labels$title)
    return(p)
  }, p = p, pane = pane, SIMPLIFY = FALSE)
}



#' @title Estimate \eqn{\hat{A}} from Tukey's \eqn{g}- or \eqn{gh}-distribution
#' 
#' @description
#' To get various \eqn{\hat{A}} from Tukey's \eqn{g}- or \eqn{gh}-distribution,
#' by minimizing the standard deviations of \eqn{\hat{g}_p}
#' 
#' @param x,probs_g see function [letterValue_g()]
#' 
#' @param A_select `'median'`, `'optim'` (default)
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous labs
#' @importFrom scales label_percent
#' @importFrom latex2exp TeX
#' @importFrom stats median.default qnorm quantile sd
#' @keywords internal
#' @export
lv_g_A <- function(x, probs_g, A_select = c('optim', 'median'), ...) {
  
  tL <- quantile(x, probs = probs_g)
  tU <- quantile(x, probs = 1 - probs_g)
  z <- qnorm(probs_g)
  
  foo <- \(A) {
    lhs <- A - tL
    uhs <- tU - A
    gs <- log(lhs / uhs) / z
    gs |> sd()
  }
  
  A_select <- match.arg(A_select)
  A_intv_ <- quantile(x, probs = c(.45, .55))
  
  A <- switch(A_select, median = {
    median.default(x)
  }, optim = {
    optimize(f = foo, interval = A_intv_)$minimum
  })
  
  gs <- log((A - tL) / (tU - A)) / z
  
  attr(A, which = 'plot') <- ggplot() +
    geom_point(mapping = aes(x = probs_g, y = gs)) +
    scale_x_continuous(name = NULL, labels = label_percent()) +
    labs(
      y = sprintf(fmt = '$\\hat{g}_p$ (sd=%.3f)', sd(gs)) |> TeX(),
      subtitle = sprintf(fmt = '$\\hat{A}_{%s}=%.3f$', A_select, A) |> TeX()
    )
  
  return(A)
  
}





#' @title Estimate \eqn{\hat{g}} from Tukey's \eqn{g}-distribution
#' 
#' @description
#' To get various \eqn{\hat{g}} from Tukey's \eqn{g}-distribution,
#' by minimizing the difference between \eqn{\hat{B}_U} and \eqn{\hat{B}_L}. 
#' 
#' @param x,probs_g see function [letterValue_g()]
#' 
#' @param A \link[base]{double} scalar, estimate \eqn{\hat{A}}
#' 
#' @param g_select `'median'`, `'optim'` (default)
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @importFrom stats optimize qnorm quantile
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_y_continuous dup_axis guide_axis labs
#' @importFrom latex2exp TeX
#' @importFrom scales label_percent
#' @keywords internal
#' @export
lv_g_g <- function(x, A, probs_g, g_select = c('optim', 'median', 'demo'), ...) {
  
  tL <- quantile(x, probs = probs_g)
  tU <- quantile(x, probs = 1 - probs_g)
  z <- qnorm(probs_g)
  gs <- log((A - tL) / (tU - A)) / z
  g_intv_ <- quantile(gs, probs = c(.4, .6))
  
  g_select <- match.arg(g_select)
  
  if (g_select %in% c('median', 'demo')) {
    g_median_ <- median.default(gs)
  }
  
  foo <- \(g) {
    BL <- .lm.fit(y = cbind(tL - A), x = cbind(expm1(g*z) / g))$coefficient
    BU <- .lm.fit(y = cbind(tU - A), x = cbind(expm1(-g*z) / g))$coefficient
    return((BL - BU)^2)
  }
  
  if (g_select %in% c('optim', 'demo')) {
    g_optim_ <- optimize(f = foo, interval = g_intv_)$minimum
  }
  
  g <- switch(g_select, median = {
    c(median = g_median_)
  }, optim = {
    c(optim = g_optim_)
  }, demo = {
    c(median = g_median_, optim = g_optim_)
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





#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous labs
#' @importFrom latex2exp TeX
#' @importFrom stats .lm.fit qnorm quantile
#' @importFrom rlang .data
lv_B_ <- function(x, A, g, probs_g) {
  q <- quantile(x, probs = c(probs_g, 1-probs_g)) - A
  z <- qnorm(probs_g)
  n <- length(probs_g)
  sq <- seq_len(n)
  r_x <- expm1(c(g*z, -g*z))/g
  B <- unname(.lm.fit(x = cbind(r_x), y = cbind(q))$coefficients)
  if (B < 0) stop('estimated B < 0')
  BL <- unname(.lm.fit(x = cbind(r_x[sq]), y = cbind(q[sq]))$coefficients)
  if (BL < 0) stop('estimated B (lower spread) < 0')
  BU <- unname(.lm.fit(x = cbind(r_x[-sq]), y = cbind(q[-sq]))$coefficients)
  if (BU < 0) stop('estimated B (upper spread) < 0')
  
  aes_d_ <- data.frame(
    x = r_x, y = q,
    color = c('lower', 'upper') |> rep(each = n), # in case `BL` and `BU` too close
    label1 = sprintf(fmt = 'B=%.3f', c(BL, BU)) |> rep(each = n),
    label2 = sprintf(fmt = 'B=%.3f', B)
  )
  
  attr(B, which = 'plot') <- ggplot() + 
    
    geom_point(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, color = .data$color), alpha = .1, show.legend = FALSE) + 
    
    geom_textsmooth(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, color = .data$color, label = .data$label1), formula = y ~ x, linewidth = .5, linetype = 0L, method = 'lm', show.legend = FALSE) +
    
    geom_textsmooth(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, label = .data$label2), color = 'grey40', formula = y ~ x, linewidth = .5, linetype = 0L, method = 'lm', show.legend = FALSE) +
    
    scale_x_continuous(
      name = '$(e^{gz_p} - 1)/g$' |> TeX(),
      sec.axis = sec_axis(
        name = 'p', 
        transform = ~ pnorm(log(.*g+1)/g),
        breaks = c(.05, .1, .25, .5, .75, .9, .95),
        labels = label_percent()
      )
    ) +
    
    labs(
      subtitle = sprintf(fmt = '$\\hat{g}_{%s} = %.3f$', names(g), g) |> TeX(),
      y = '$t_p$' |> TeX(), 
      caption = 'Constrained at h = 0'
    )
  
  return(B)
}




