# FUNCTION: evaluate the probablity of the mixture distribution
# inputs: y is a vector of observations (of length N)
# a_mcmc, b_mcmc, p_mcmc are vectors of length L
# outputs: a vector of probablity values that has the same dimension as y
# first outputs a NxL matrix for the probablity evaluated under each of the L kernels
# then take weighted sum to get the overall probability
eval_mixture_grouped <- function(y, p_mcmc, a_mcmc, b_mcmc) {
  
  tmat <- matrix(unlist(sapply(
    y,
    FUN = "dnsbeta",
    shape1 = a_mcmc,
    shape2 = b_mcmc,
    min = 0, max = T_kernel,
  )), nrow = length(y), byrow = TRUE)
  tmat[tmat == Inf] <- 0
  temp <- sweep(tmat, MARGIN = 2, p_mcmc, `*`)
  return(apply(temp, MARGIN = 1, FUN = "sum"))
}
