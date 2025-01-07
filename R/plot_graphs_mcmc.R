# FUNCTION: plot result graph
# FUNCTION: GENERATE GRAPHS -- GENERATE SAMPLES
plot_graphs_mcmc <- function(p_mcmc, as_mcmc, bs_mcmc, p0_mcmc, a0_mcmc, b0_mcmc, epsilon_mcmc) {
  PLOT_START <- 0.01
  x_plot <- seq(PLOT_START, 1-PLOT_START, by = PLOT_START)
  n_xplot <- length(x_plot)
  B_iter <- dim(p0_mcmc)[1]
  B_burn <- floor(B_iter/2) - 1
  
  ## calculate true density on a grid 
  # 1. common component part
  pos_p0s_true <- NA
  for (j in 1:length(x_plot)) {
    pos_p0s_true[j] <- dnsbeta(x_plot[j], a0s, b0s, 0, T_kernel)
  }
  # 2. heterogenous components part 
  pos_ps_true <- list()
  pos_ps_true_mixed <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      temp <- NA
      for (j in 1:length(x_plot)) {
        temp[j] <- ps[[id_i]][[id_j]][1]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][1], bs[[id_i]][[id_j]][1], 0, T_kernel) +
          ps[[id_i]][[id_j]][2]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][2], bs[[id_i]][[id_j]][2], 0, T_kernel)
      }
      pos_ps_true[[(id_i - 1)*K + id_j]] <- temp
      pos_ps_true_mixed[[(id_i - 1)*K + id_j]] <- pos_p0s_true*(epsilon_true) + pos_ps_true[[(id_i - 1)*K + id_j]]*(1 - epsilon_true)
    }
  }
  
  ## calculate posterior results 
  pos_p0s <- matrix(0, nrow = B_iter - B_burn + 1, ncol = length(x_plot))
  for (j in 1:length(x_plot)) {
    pos_p0s[,j] <- eval_mixture(x_plot[j], 
                                p0_mcmc[B_burn:B_iter, ], 
                                a0_mcmc[B_burn:B_iter, ],
                                b0_mcmc[B_burn:B_iter, ])
  }
  
  # idio part
  pos_ps <- list()
  # mixed part
  pos_ps_mixed <- list()
  
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos_ps[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter - B_burn + 1, ncol = length(x_plot))
      temp <- NA
      for (j in 1:length(x_plot)) {
        pos_ps[[(id_i - 1)*K + id_j]][,j] <- eval_mixture(x_plot[j], p_mcmc[[(id_i - 1)*K + id_j]][B_burn:B_iter, ],
                                                          as_mcmc[[(id_i - 1)*K + id_j]][B_burn:B_iter, ],
                                                          bs_mcmc[[(id_i - 1)*K + id_j]][B_burn:B_iter, ])
      }
      pos_ps_mixed[[(id_i - 1)*K + id_j]] <- 
        pos_p0s * matrix(rep(epsilon_mcmc[B_burn:B_iter], n_xplot), ncol = n_xplot) + 
        pos_ps[[(id_i - 1)*K + id_j]] * matrix(rep((1-epsilon_mcmc[B_burn:B_iter]), n_xplot), ncol = n_xplot)
    }
  }
  colors <- c("#F76D5E", "#FFFFBF", "#72D8FF", "purple", "grey")
  plist <- list()
  for (i in 1:4) {
    df_plot_p <- data.frame(
      x = x_plot,
      mcmc = pos_ps_mixed[[i]] |> colMedians(),
      true = pos_p0s_true*(epsilon_true) + pos_ps_true[[i]]*(1 - epsilon_true)
    ) %>% pivot_longer(cols = 2:3)
    plist[[i]] <- ggplot(df_plot_p, aes(x = x, y = value, color = name)) + geom_line() +
      geom_area(aes(fill = name, group = name),
                alpha = 0.5, position = 'identity') +
      scale_fill_manual(values = c(colors[1], colors[3])) +
      ggtitle(paste0("Mixed: ", floor((i-1)/2)+1, " -> " ,((i-1) %% 2) + 1)) +
      ylab("density")
  }
  for (i in 1:4) {
    df_plot_p <- data.frame(
      x = x_plot,
      mcmc = pos_ps[[i]] |> colMedians(),
      true = pos_ps_true[[i]]
    ) %>% pivot_longer(cols = 2:3)
    plist[[i+4]] <- ggplot(df_plot_p, aes(x = x, y = value, color = name)) + geom_line() +
      geom_area(aes(fill = name, group = name),
                alpha = 0.5, position = 'identity') +
      scale_fill_manual(values = c(colors[1], colors[3])) +
      ggtitle(paste0("Idio: ",floor((i-1)/2)+1, " -> " ,((i-1) %% 2) + 1)) +
      ylab("density")
  }
  df_plot_p <- data.frame(
    x = x_plot,
    mcmc = pos_p0s |> colMedians(),
    true = pos_p0s_true
  ) %>% pivot_longer(cols = 2:3)
  plist[[9]] <- ggplot(df_plot_p, aes(x = x, y = value, color = name)) + geom_line() +
    geom_area(aes(fill = name, group = name),
              alpha = 0.5, position = 'identity') +
    scale_fill_manual(values = c(colors[1], colors[3])) +
    ggtitle(paste0("Common: ",floor((i-1)/2)+1, " -> " ,((i-1) %% 2) + 1)) +
    ylab("density")
  return(plist)
}
