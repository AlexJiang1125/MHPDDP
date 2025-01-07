calc_rmse_vi_mean <- function(alpha0_vi, alphas_vi, cas_vi, cbs_vi, das_vi, dbs_vi,
                              ca0_vi, cb0_vi, da0_vi, db0_vi, eps_vi, type = "RANDOM") {
  res_list <- get_DPpars(EXAMPLE_ID = EXAMPLE_ID, samepriors = samepriors, epsilon_true = eps_true)
  list2env(res_list, envir = globalenv())
  res_list <- get_HPpars(HP_ID = HP_ID)
  list2env(res_list, envir = globalenv())
  
  PLOT_START <- 0.01
  x_plot <- seq(PLOT_START, 1-PLOT_START, by = PLOT_START)
  n_xplot <- length(x_plot)
  # generate samples from the variational distribution
  p0m <- alpha0_vi/sum(alpha0_vi)
  a0m <- ca0_vi/da0_vi
  b0m <- cb0_vi/db0_vi
  psm <- asm <- bsm <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      psm[[(id_i - 1)*K + id_j]] <- alphas_vi[[(id_i - 1)*K + id_j]]/sum(alphas_vi[[(id_i - 1)*K + id_j]])
      asm[[(id_i - 1)*K + id_j]] <- cas_vi[[(id_i - 1)*K + id_j]]/das_vi[[(id_i - 1)*K + id_j]]
      bsm[[(id_i - 1)*K + id_j]] <- cbs_vi[[(id_i - 1)*K + id_j]]/dbs_vi[[(id_i - 1)*K + id_j]]
    }
  }
  if (type == "RANDOM") {
    epm <- eps_vi[1]/sum(eps_vi)
  } else if (type == "COMMON") {
    epm <- 1
  } else {
    epm <- 0
  }
  pos_p0s <- NA
  pos_p0s_true <- NA
  for (j in 1:length(x_plot)) {
    pos_p0s[j] <- sum(dbeta(x_plot[j], a0m, b0m)*p0m)
    pos_p0s_true[j] <- dbeta(x_plot[j], a0s, b0s)
  }
  pos_p0s_median <- pos_p0s
  # normalize: make sure that integrates to one
  pos_p0s_median <- pos_p0s_median/sum(PLOT_START*pos_p0s_median)
  colors <- c("#F76D5E", "#FFFFBF", "#72D8FF", "purple", "grey")
  pos_ps <- list()
  pos_ps_median <- list()
  pos_ps_true <- list()
  pos_ps_mixed <- list()
  pos_ps_mixed_median <- list()
  RMISE <- NA
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos_ps[[(id_i - 1)*K + id_j]] <- NA
      temp <- NA
      for (j in 1:length(x_plot)) {
        pos_ps[[(id_i - 1)*K + id_j]][j] <- sum(dbeta(x_plot[j],
                                                      asm[[(id_i - 1)*K + id_j]],
                                                      bsm[[(id_i - 1)*K + id_j]])*psm[[(id_i - 1)*K + id_j]])
        temp[j] <- ps[[id_i]][[id_j]][1]*dbeta(x_plot[j], as[[id_i]][[id_j]][1], bs[[id_i]][[id_j]][1]) +
          ps[[id_i]][[id_j]][2]*dbeta(x_plot[j], as[[id_i]][[id_j]][2], bs[[id_i]][[id_j]][2])
      }
      pos_ps_true[[(id_i - 1)*K + id_j]] <- temp
      pos_ps_median[[(id_i - 1)*K + id_j]] <- pos_ps[[(id_i - 1)*K + id_j]]
      pos_ps_median[[(id_i - 1)*K + id_j]] <- pos_ps_median[[(id_i - 1)*K + id_j]]/sum(PLOT_START*pos_ps_median[[(id_i - 1)*K + id_j]])
      pos_ps_mixed[[(id_i - 1)*K + id_j]] <-
        pos_p0s * epm +
        pos_ps[[(id_i - 1)*K + id_j]] * (1-epm)
      pos_ps_mixed_median[[(id_i - 1)*K + id_j]] <- pos_ps_mixed[[(id_i - 1)*K + id_j]]
      pos_ps_mixed_median[[(id_i - 1)*K + id_j]] <- pos_ps_mixed_median[[(id_i - 1)*K + id_j]]/sum(PLOT_START*pos_ps_mixed_median[[(id_i - 1)*K + id_j]])
      pos_ps_true_mixed <- pos_p0s_true*(epsilon_true) + pos_ps_true[[(id_i - 1)*K + id_j]]*(1 - epsilon_true)
      
      RMISE[(id_i - 1)*K + id_j] <- sqrt(sum(PLOT_START*(pos_ps_mixed_median[[(id_i - 1)*K + id_j]] - pos_ps_true_mixed)^2))
    }
  }
  return(RMISE)
}
