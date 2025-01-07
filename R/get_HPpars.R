get_HPpars <- function(HP_ID, flipped = FALSE) {

  K <- 2
  e_alpha <- 2
  f_alpha <- 4

  if (HP_ID == 4) {
    T_all <- 150000
    a_mu <- 1
    b_mu <- 15
  } else {
    T_all <- 15000
    a_mu <- 1
    b_mu <- 100
  }

  if (HP_ID == 1) {
    alpha_true <- matrix(c(0.6, 0.15, 0.3, 0.6), nrow = 2)
    mu_true<- c(0.05, 0.1)
  } else if (HP_ID == 2) {
    alpha_true <- matrix(c(0.6, 0.3, 0.3, 0.6), nrow = 2)
    mu_true<- c(0.005, 0.005)
  } else if (HP_ID == 3) {
    alpha_true <- matrix(c(0.4, 0.2, 0.1, 0.4), nrow = 2)
    mu_true<- c(0.02, 0.02)
  } else if (HP_ID == 4) {
    alpha_true <- matrix(c(0.7, 0.05, 0.3, 0.6), nrow = 2)
    mu_true<- c(0.01, 0.01)
  }
  return(
    res = list(
      alpha_true = alpha_true,
      mu_true = mu_true,
      T_all = T_all, K = K,
      a_mu = a_mu, b_mu = b_mu,
      e_alpha = e_alpha, f_alpha = f_alpha
    )
  )
}


