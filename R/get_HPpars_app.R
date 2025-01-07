get_HPpars_app <- function(HP_ID, filename, flipped = FALSE) {

  if (filename == "app_buy_orders_byvol.csv" || filename == "app_sell_orders_byvol.csv") {
    K = 4
  } else if (filename == "app_orders_all.csv") {
    K = 2
  } else {
    K = 4
  }
  e_alpha <- 2
  f_alpha <- 4

  T_all = 2340
  a_mu <- 1
  b_mu <- 15

  alpha_true <- matrix(c(0.6, 0.15, 0.3, 0.6), nrow = 2)
  mu_true<- c(0.05, 0.1)

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


