# FUNCTION: evaluate the probablity of the mixture distribution
eval_mixture <- function(y, p_mcmc, a_mcmc, b_mcmc, T_kernel = 1, type = "d") {
  if (type == "d") {
    return(apply(matrix(dnsbeta(y, a_mcmc, b_mcmc, 0, T_kernel), nrow = dim(a_mcmc)[1],
                        ncol = dim(a_mcmc)[2])*p_mcmc, FUN = "sum", MARGIN = 1))
  } else {
    return(apply(matrix(pnsbeta(y, a_mcmc, b_mcmc, 0, T_kernel), nrow = dim(a_mcmc)[1],
                        ncol = dim(a_mcmc)[2])*p_mcmc, FUN = "sum", MARGIN = 1))
  }
}

my_dnsbeta <- function(x, a_mean, b_mean, min = 0, max = T_kernel, log = FALSE) {
  a_shape <- a_mean/(1-b_mean)
  b_shape <- a_mean/b_mean
  return(dnsbeta(x, a_shape, b_shape, min = min, max = max, log = log))
}

my_pnsbeta <- function(q, a_mean, b_mean, min = 0, max = T_kernel) {
  a_shape <- a_mean/(1-b_mean)
  b_shape <- a_mean/b_mean
  return(pnsbeta(q = q, shape1 = a_shape, shape2 = b_shape, min = min, max = max))
}

eval_mixture_meanspec <- function(y, p_mcmc, a_mcmc, b_mcmc, type = "d") {
  if (type == "d") {
    return(apply(matrix(my_dnsbeta(y, a_mcmc, b_mcmc, 0, T_kernel), nrow = dim(b_mcmc)[1],
                        ncol = dim(a_mcmc)[2])*p_mcmc, FUN = "sum", MARGIN = 1))
  } else {
    return(apply(matrix(my_pnsbeta(y, a_mcmc, b_mcmc, 0, T_kernel), nrow = dim(b_mcmc)[1],
                        ncol = dim(a_mcmc)[2])*p_mcmc, FUN = "sum", MARGIN = 1))
  }
}
