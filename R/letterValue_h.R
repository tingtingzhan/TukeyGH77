

#' @rdname letterValue
#' @import patchwork
#' @importFrom geomtextpath geom_textsmooth
#' @export
letterValue_h <- function(
    x,
    probs_h = seq.int(from = .1, to = .3, by = .005),
    A_select = c('median', 'B.optim', 'h.optim', 'demo'),
    #do.plot = FALSE,
    ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  
  probs_h <- probs_h |> check_letterVal_(x = x)
  
  z <- qnorm(probs_h)
  r_x <- z^2 / 2 # for final plot!
  r_X <- cbind(1, r_x) # for ?stats::.lm.fit
  t. <- quantile(x, probs = c(1 - probs_h, probs_h))
  sq <- seq_along(probs_h)
  r_y <- log((t.[sq] - t.[-sq]) / (-2*z)) # for final plot!
  cf <- unname(.lm.fit(x = r_X, y = cbind(r_y))$coefficients)
  B <- exp(cf[1L])
  h <- cf[2L]
  Bh_layer <- geom_textsmooth(mapping = aes(x = r_x, y = r_y), label = sprintf(fmt = 'B=%.3f, h=%.3f', B, h), color = 'grey40', formula = y ~ x, linewidth = .5, linetype = 2L, method = 'lm', show.legend = FALSE)
  # has nothing to do \hat{A} !!!
  
  A_select <- match.arg(A_select)
  
  if (A_select == 'demo') {# plot only!!
    p <- mapply(FUN = \(s, pane) {
      p <- lv_h_A(x = x, probs_h = probs_h, A_select = s) |>
        attr(which = 'plot', exact = TRUE)
      p$labels$title <- sprintf(fmt = '(%s)', pane) |> 
        paste(p$labels$title)
      p + Bh_layer
    }, s = c('median', 'B.optim', 'h.optim'), pane = LETTERS[seq_len(3L)], SIMPLIFY = FALSE) |>
      Reduce(f = `+`)
    return(p)
  }
  
  A <- lv_h_A(x = x, probs_h = probs_h, A_select = A_select)
  ret <- c(A = A, B = B, g = 0, h = pmax(0, h))
  #if (do.plot) attr(ret, which = 'plot') <- attr(A, which = 'plot', exact = TRUE) + Bh_layer
  return(ret)
  
}


#' @title Estimate \eqn{\hat{A}} from Tukey's \eqn{h}-distribution
#' 
#' @description
#' To get various \eqn{\hat{A}} from Tukey's \eqn{h}-distribution,
#' by optimizing letter-value based estimates on lower- and upper- half spread.
#' 
#' @param x,probs_h see function [letterValue_h()]
#' 
#' @param A_select `'median'` (default), `'B.optim'` or `'h.optim'`
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @importFrom geomtextpath geom_textabline geom_textsmooth
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous sec_axis labs
#' @importFrom stats .lm.fit median.default qnorm quantile
#' @importFrom scales label_percent
#' @importFrom rlang .data
#' @keywords internal
#' @export
lv_h_A <- function(x, probs_h, A_select = c('median', 'B.optim', 'h.optim'), ...) {

  z <- qnorm(probs_h)
  r_x <- z^2 / 2
  r_X <- cbind(1, r_x) # for ?stats::.lm.fit
  tL <- quantile(x, probs = probs_h)
  tU <- quantile(x, probs = 1 - probs_h)
  n <- length(probs_h)
  sq <- seq_len(n)
  
  foo <- \(A) {
    lhs <- A - tL
    uhs <- tU - A
    r_y <- log(-c(lhs, uhs) / z)
    dim(r_y) <- c(2*n, 1L)
    cfL <- unname(.lm.fit(x = r_X, y = r_y[sq, , drop = FALSE])$coefficients)
    cfU <- unname(.lm.fit(x = r_X, y = r_y[-sq, , drop = FALSE])$coefficients)
    ret <- (cfL - cfU)^2
    names(ret) <- c('errB', 'errh')
    return(ret)
  }
  
  A_select <- match.arg(A_select)
  A_intv_ <- quantile(x, probs = c(.45, .55))
  
  A <- switch(A_select, median = {
    median.default(x)
  }, B.optim = {
    optimize(f = \(A) foo(A)['errB'], interval = A_intv_)$minimum
  }, h.optim = {
    optimize(f = \(A) foo(A)['errh'], interval = A_intv_)$minimum
  })
  
  r_yL <- log(-(A - tL)/z)
  r_yU <- log(-(tU - A)/z)
  cfL <- unname(.lm.fit(x = r_X, y = cbind(r_yL))$coefficients)
  cfU <- unname(.lm.fit(x = r_X, y = cbind(r_yU))$coefficients)
  
  aes_d_ <- data.frame(
    x = r_x,
    y = c(r_yL, r_yU), 
    label = sprintf(fmt = '%s: B=%.3f, h=%.3f', c('Lower', 'Upper'), exp(c(cfL[1L], cfU[1L])), c(cfL[2L], cfU[2L])) |> rep(each = n)
  )

  attr(A, which = 'plot') <- ggplot() + 
    
    geom_point(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, color = .data$label), alpha = .2, show.legend = FALSE) + 
    
    geom_textsmooth(data = aes_d_, mapping = aes(x = .data$x, y = .data$y, label = .data$label, color = .data$label), formula = y ~ x, linewidth = .5, linetype = 2L, method = 'lm', show.legend = FALSE) +
    
    scale_x_continuous(
      name = '$z_p^2/2$' |> TeX(),
      sec.axis = sec_axis(
        name = 'p', 
        transform = ~ pnorm(sqrt(. * 2), lower.tail = FALSE),
        labels = label_percent()
      )
    ) +
    
    labs(
      subtitle = sprintf(fmt = '$\\hat{A}_{%s}=%.3f$', A_select, A) |> TeX(),
      y = NULL,
      caption = 'Constrained at g = 0'
    )
  
  return(A)
  
}



