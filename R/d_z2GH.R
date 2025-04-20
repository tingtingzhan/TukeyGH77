
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





