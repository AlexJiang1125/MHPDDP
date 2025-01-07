# FUNCTION: weight average
weight_average <- function(x_old, x_new, rho) {
  return(x_new*rho + x_old*(1-rho))
}
