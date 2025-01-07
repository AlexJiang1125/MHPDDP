# load the parameters for Dirichlet processes
get_DPpars_app <- function(EXAMPLE_ID, T_kernel = 1) {
  alpha_DP <- 1 # DP concentration parameter
  #T_kernel <- 1
  c_a <- c_a0 <- 0.5
  d_a <- d_a0 <- 1
  
  c_b <- c_b0 <- 2
  d_b <- d_b0 <- 1
  
  return(
    res = list(
      c_a = c_a, c_b = c_b,
      c_a0 = c_a0, c_b0 = c_b0,
      d_a = d_a, d_b = d_b,
      d_a0 = d_a0, d_b0 = d_b0,
      alpha_DP = alpha_DP,
      T_kernel = T_kernel
    )
  )
}

