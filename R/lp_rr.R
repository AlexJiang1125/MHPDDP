# FUNCTION: log-prior part using prior specifications from Robert and Rosseau (2002)
lp_rr <- function(s,m) {
  sq_sum <- sum((c(s,m) - c(2,0.5))^2)
  part_1 <- log(1 - exp(-delta*(sq_sum)))
  part_2 <- -rho/s^2/m/(1-m) - kappa*s^2/2
  return(part_1 + part_2)
}
