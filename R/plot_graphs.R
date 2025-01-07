# FUNCTION: plot result graph
# FUNCTION: GENERATE GRAPHS -- GENERATE SAMPLES
plot_graphs <- function(alpha0_vi, alphas_vi, cas_vi, cbs_vi, das_vi, dbs_vi,
                        ca0_vi, cb0_vi, da0_vi, db0_vi, eps_vi) {
  B_burn <- floor(B/2) + 1
  PLOT_START <- 0.01
  x_plot <- seq(PLOT_START, 1-PLOT_START, by = PLOT_START)
  n_xplot <- length(x_plot)
  # generate samples from the variational distribution
  p0_mcmc <- rdirichlet(n = B-B_burn+1, alpha = alpha0_vi)
  a0_mcmc <- matrix(0, nrow = B-B_burn+1, ncol = L)
  b0_mcmc <- matrix(0, nrow = B-B_burn+1, ncol = L)
  for (l in 1:L) {
    a0_mcmc[, l] <- rgamma(n = B-B_burn+1, ca0_vi[l], da0_vi[l])
    b0_mcmc[, l] <- rgamma(n = B-B_burn+1, cb0_vi[l], db0_vi[l])
  }
  p_mcmc <- as_mcmc <- bs_mcmc <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      p_mcmc[[(id_i - 1)*K + id_j]] <- rdirichlet(n = B-B_burn+1, alpha = alpha0_vi)
      as_mcmc[[(id_i - 1)*K + id_j]] <- rgamma(n = B-B_burn+1, cas_vi[[(id_i - 1)*K + id_j]][l], das_vi[[(id_i - 1)*K + id_j]][l])
      bs_mcmc[[(id_i - 1)*K + id_j]] <- rgamma(n = B-B_burn+1, cbs_vi[[(id_i - 1)*K + id_j]][l], dbs_vi[[(id_i - 1)*K + id_j]][l])
    }
  }
  epsilon_mcmc <- rnsbeta(n = B-B_burn+1, eps_vi[1], eps_vi[2], 0, T_kernel)
  pos_p0s <- matrix(0, nrow = B - B_burn + 1, ncol = length(x_plot))
  pos_p0s_true <- NA
  for (j in 1:length(x_plot)) {
    pos_p0s[,j] <- eval_mixture(x_plot[j], p0_mcmc, a0_mcmc, b0_mcmc)
    pos_p0s_true[j] <- dnsbeta(x_plot[j], a0s, b0s, 0, T_kernel)
  }
  pos_p0s_median <- apply(pos_p0s, FUN = "median", 2)
  # normalize: make sure that integrates to one
  pos_p0s_median <- pos_p0s_median/sum(PLOT_START*pos_p0s_median)
  colors <- c("#F76D5E", "#FFFFBF", "#72D8FF", "purple", "grey")
  pos_ps <- list()
  pos_ps_median <- list()
  pos_ps_true <- list()
  pos_ps_mixed <- list()
  pos_ps_mixed_median <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos_ps[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B - B_burn + 1, ncol = length(x_plot))
      temp <- NA
      for (j in 1:length(x_plot)) {
        pos_ps[[(id_i - 1)*K + id_j]][,j] <- eval_mixture(x_plot[j], p_mcmc[[(id_i - 1)*K + id_j]],
                                                          as_mcmc[[(id_i - 1)*K + id_j]],
                                                          bs_mcmc[[(id_i - 1)*K + id_j]])
        temp[j] <- ps[[id_i]][[id_j]][1]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][1], bs[[id_i]][[id_j]][1], 0, T_kernel) +
          ps[[id_i]][[id_j]][2]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][2], bs[[id_i]][[id_j]][2], 0, T_kernel)
      }
      pos_ps_true[[(id_i - 1)*K + id_j]] <- temp
      pos_ps_median[[(id_i - 1)*K + id_j]] <- apply(pos_ps[[(id_i - 1)*K + id_j]], FUN = "median", 2)
      pos_ps_median[[(id_i - 1)*K + id_j]] <- pos_ps_median[[(id_i - 1)*K + id_j]]/sum(PLOT_START*pos_ps_median[[(id_i - 1)*K + id_j]])
      pos_ps_mixed[[(id_i - 1)*K + id_j]] <-
        pos_p0s * matrix(rep(epsilon_mcmc, n_xplot), ncol = n_xplot) +
        pos_ps[[(id_i - 1)*K + id_j]] * matrix(rep((1-epsilon_mcmc), n_xplot), ncol = n_xplot)
      pos_ps_mixed_median[[(id_i - 1)*K + id_j]] <- apply(pos_ps_mixed[[(id_i - 1)*K + id_j]], FUN = "median", 2)
      pos_ps_mixed_median[[(id_i - 1)*K + id_j]] <- pos_ps_mixed_median[[(id_i - 1)*K + id_j]]/sum(PLOT_START*pos_ps_mixed_median[[(id_i - 1)*K + id_j]])
    }
  }
  plist <- list()
  for (i in 1:4) {
    df_plot_p <- data.frame(
      x = x_plot,
      mcmc = pos_ps_mixed_median[[i]],
      true = pos_p0s_true*(epsilon_true) + pos_ps_true[[i]]*(1 - epsilon_true)
    ) %>% pivot_longer(cols = 2:3)
    plist[[i]] <- ggplot(df_plot_p, aes(x = x, y = value, color = name)) + geom_line() +
      geom_area(aes(fill = name, group = name),
                alpha = 0.5, position = 'identity') +
      scale_fill_manual(values = c(colors[1], colors[3])) +
      ggtitle(paste0(floor((i-1)/2)+1, " -> " ,((i-1) %% 2) + 1)) +
      ylab("density")
  }
  for (i in 1:4) {
    df_plot_p <- data.frame(
      x = x_plot,
      mcmc = pos_ps_median[[i]],
      true = pos_ps_true[[i]]
    ) %>% pivot_longer(cols = 2:3)
    plist[[i+4]] <- ggplot(df_plot_p, aes(x = x, y = value, color = name)) + geom_line() +
      geom_area(aes(fill = name, group = name),
                alpha = 0.5, position = 'identity') +
      scale_fill_manual(values = c(colors[1], colors[3])) +
      ggtitle(paste0(floor((i-1)/2)+1, " -> " ,((i-1) %% 2) + 1)) +
      ylab("density")
  }
  df_plot_p <- data.frame(
    x = x_plot,
    mcmc = pos_p0s_median,
    true = pos_p0s_true
  ) %>% pivot_longer(cols = 2:3)
  plist[[9]] <- ggplot(df_plot_p, aes(x = x, y = value, color = name)) + geom_line() +
    geom_area(aes(fill = name, group = name),
              alpha = 0.5, position = 'identity') +
    scale_fill_manual(values = c(colors[1], colors[3])) +
    ggtitle(paste0(floor((i-1)/2)+1, " -> " ,((i-1) %% 2) + 1)) +
    ylab("density")
  return(plist)
}
