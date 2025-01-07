calc_rmse_mcmc <- function(p_mcmc, as_mcmc, bs_mcmc, p0_mcmc, a0_mcmc, b0_mcmc, epsilon_mcmc, B_iter, fig_dir) {
  library(EnvStats)

  interval_score <- function(true_values, lower, upper, interval_range, weigh = TRUE,
                             separate_results = FALSE)
  {
    present <- c(methods::hasArg("true_values"), methods::hasArg("lower"),
                 methods::hasArg("upper"), methods::hasArg("interval_range"))
    if (!all(present)) {
      stop("need all arguments 'true_values', 'lower', 'upper' and 'interval_range' in function 'interval_score()'")
    }
    alpha <- (100 - interval_range)/100
    dispersion <- (upper - lower)
    overprediction <- 2/alpha * (lower - true_values) * (true_values <
                                                           lower)
    underprediction <- 2/alpha * (true_values - upper) * (true_values >
                                                            upper)
    if (weigh) {
      dispersion <- dispersion * alpha/2
      underprediction <- underprediction * alpha/2
      overprediction <- overprediction * alpha/2
    }
    score <- dispersion + underprediction + overprediction
    if (separate_results) {
      return(list(interval_score = score, dispersion = dispersion,
                  underprediction = underprediction, overprediction = overprediction))
    }
    else {
      return(score)
    }
  }

  #library(scoringutils)
  res_list <- get_DPpars(EXAMPLE_ID = EXAMPLE_ID, samepriors = samepriors, epsilon_true = eps_true)
  list2env(res_list, envir = globalenv())
  res_list <- get_HPpars(HP_ID = HP_ID)
  list2env(res_list, envir = globalenv())

  x_plot <- seq(PLOT_START, T_kernel-PLOT_START, by = PLOT_START)
  n_xplot <- length(x_plot)
  B_iter <- length(epsilon_mcmc)
  B_burn <- floor(B_iter/2) - 1

  if (EXAMPLE_ID == 8) {
    T_kernel = 5
  } else {
    T_kernel = 1
  }

  ## calculate true density on a grid
  # 1. common component part
  pos_p0s_true <- NA
  if (EXAMPLE_ID == 5) {
    for (j in 1:length(x_plot)) {
      pos_p0s_true[j] <- eps_tri0*dtri(x_plot[j], min = tri_a0[1], max = tri_c0[1], mode = tri_b0[1]) +
        (1-eps_tri0)*dtri(x_plot[j], min = tri_a0[2], max = tri_c0[2], mode = tri_b0[2])
    }
  } else if (EXAMPLE_ID %in% c(6,8)) {
    for (j in 1:length(x_plot)) {
      pos_p0s_true[j] <- dexp(x_plot[j], exp_0)
    }
  } else if (EXAMPLE_ID == 7) {
    length_p0mixtures <- length(p0s)
    for (j in 1:length(x_plot)) {
      for (id_l in 1:length_p0mixtures) {
        if (id_l == 1) {
          pos_p0s_true[j] <- p0s[1]*dnsbeta(x_plot[j], a0s[1], b0s[1], 0, T_kernel)
        } else {
          pos_p0s_true[j] <- pos_p0s_true[j] + p0s[id_l]*dnsbeta(x_plot[j], a0s[id_l], b0s[id_l], 0, T_kernel)
        }
      }
    }
  } else {
    for (j in 1:length(x_plot)) {
      pos_p0s_true[j] <- dnsbeta(x_plot[j], a0s, b0s,0,T_kernel)
    }
  }
  # 2. heterogenous components part
  pos_ps_true <- list()
  pos_ps_true_mixed <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos <- (id_i - 1)*K + id_j
      temp <- NA
      for (j in 1:length(x_plot)) {
        length_pmixtures <- length(ps[[id_i]][[id_j]])
        if (EXAMPLE_ID == 5) {
          if (pos == 4) {
            temp[j] <- dnsbeta(x_plot[j], 1, 1, 0, T_kernel)
          } else {
            temp[j] <- eps_tri*dtri(x_plot[j], min = tri_a[[pos]][1], max = tri_c[[pos]][1], mode = tri_b[[pos]][1]) +
              (1-eps_tri)*dtri(x_plot[j], min = tri_a[[pos]][2], max = tri_c[[pos]][2], mode = tri_b[[pos]][2])
          }
        } else if (EXAMPLE_ID %in% c(6,8)) {
          temp[j] <- dexp(x_plot[j], exp_s[[id_i]][[id_j]])
        } else {
          # for examples 4 (single mixture) and 7 (DP mixture)
          for (id_l in 1:length_pmixtures) {
            if (id_l == 1) {
              temp[j] <- ps[[id_i]][[id_j]][1]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][1], bs[[id_i]][[id_j]][1], 0, T_kernel)
            } else {
              temp[j] <- temp[j] + ps[[id_i]][[id_j]][id_l]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][id_l], bs[[id_i]][[id_j]][id_l], 0, T_kernel)
            }
          }
        }
      }
      pos_ps_true[[(id_i - 1)*K + id_j]] <- temp
      pos_ps_true_mixed[[(id_i - 1)*K + id_j]] <- pos_p0s_true*(eps_true) + pos_ps_true[[(id_i - 1)*K + id_j]]*(1 - eps_true)
    }
  }

  ## calculate posterior results
  pos_p0s <- matrix(0, nrow = B_iter - B_burn + 1, ncol = length(x_plot))
  for (j in 1:length(x_plot)) {
    pos_p0s[,j] <- eval_mixture(x_plot[j],
                                p0_mcmc[B_burn:B_iter, ],
                                a0_mcmc[B_burn:B_iter, ],
                                b0_mcmc[B_burn:B_iter, ], T_kernel = T_kernel)
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
                                                          bs_mcmc[[(id_i - 1)*K + id_j]][B_burn:B_iter, ], T_kernel = T_kernel)
      }
      pos_ps_mixed[[(id_i - 1)*K + id_j]] <-
        pos_p0s * matrix(rep(epsilon_mcmc[B_burn:B_iter], n_xplot), ncol = n_xplot) +
        pos_ps[[(id_i - 1)*K + id_j]] * matrix(rep((1-epsilon_mcmc[B_burn:B_iter]), n_xplot), ncol = n_xplot)
    }
  }
  # calculate
  rmse_list <- list()



  for (id_i in 1:K) {
    for (id_j in 1:K) {
      rmse_list[[(id_i - 1)*K + id_j]] <- PLOT_START*(sweep(x = pos_ps_mixed[[(id_i - 1)*K + id_j]],
                                                            MARGIN = 2,
                                                            STATS = pos_ps_true_mixed[[(id_i - 1)*K + id_j]]))^2
    }
  }
  rmse_agg <- lapply(rmse_list, FUN = "rowSums")
  RMISE <- NA
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      RMISE[(id_i - 1)*K + id_j] <- mean(sqrt(rmse_agg[[(id_i - 1)*K + id_j]]))
    }
  }

  ### pointwise credible interval
  pws_cvg <- NA
  cw <- NA
  is <- NA

  for (id_i in 1:K) {
    for (id_j in 1:K) {
      # check upper bound
      ub_check <- sweep(x = pos_ps_mixed[[(id_i - 1)*K + id_j]],
                        MARGIN = 2,
                        STATS = pos_ps_true_mixed[[(id_i - 1)*K + id_j]]) |>
        apply(FUN = "quantile", probs = 0.95, MARGIN = 2) > 0

      lb_check <- sweep(x = pos_ps_mixed[[(id_i - 1)*K + id_j]],
                        MARGIN = 2,
                        STATS = pos_ps_true_mixed[[(id_i - 1)*K + id_j]]) |>
        apply(FUN = "quantile", probs = 0.05, MARGIN = 2) < 0
      pws_cvg[(id_i - 1)*K + id_j] <- mean(ub_check & lb_check)

      cw[(id_i - 1)*K + id_j] <- mean(apply(pos_ps_mixed[[(id_i - 1)*K + id_j]], FUN = "quantile", probs = 0.95, MARGIN = 2) -
                                        apply(pos_ps_mixed[[(id_i - 1)*K + id_j]], FUN = "quantile", probs = 0.05, MARGIN = 2))

      is[(id_i - 1)*K + id_j] <- mean(interval_score(
        true_values = pos_ps_true_mixed[[(id_i - 1)*K + id_j]],
        lower = apply(pos_ps_mixed[[(id_i - 1)*K + id_j]], FUN = "quantile", probs = 0.05, MARGIN = 2),
        upper = apply(pos_ps_mixed[[(id_i - 1)*K + id_j]], FUN = "quantile", probs = 0.95, MARGIN = 2),
        interval_range = 90,
        weigh = TRUE,
        separate_results = FALSE
      ))
    }
  }

  ggplot2::theme_set(theme_bw())
  xgrid <- seq(from = 0, to = 1, length.out = 101)

  p_common <- data.frame(
    x = x_plot,
    y_true = pos_p0s_true,
    y_mcmc = colMeans(pos_p0s)
  ) %>% pivot_longer(cols = 2:3) %>% dplyr::rename(models = "name", density = "value") %>% ggplot() +
    geom_area(aes(x = x, y = density, fill = models), alpha = 0.6, position = "identity") +
    ggtitle(paste0("common")) + theme(legend.position="bottom")

  ps <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      ps[[(id_i-1)*K + id_j]] <- data.frame(
        x = x_plot,
        y_true = pos_ps_true[[(id_i-1)*K + id_j]],
        y_mcmc = colMeans(pos_ps[[(id_i-1)*K + id_j]])
      ) %>% pivot_longer(cols = 2:3) %>% dplyr::rename(models = "name", density = "value") %>% ggplot() +
        geom_area(aes(x = x, y = density, fill = models), alpha = 0.6, position = "identity") +
        ggtitle(paste0(id_i, "->", id_j)) + theme(legend.position="bottom")
    }
  }

  pms <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pms[[(id_i-1)*K + id_j]] <- data.frame(
        x = x_plot,
        y_true = pos_ps_true_mixed[[(id_i-1)*K + id_j]],
        y_mcmc = colMeans(pos_ps_mixed[[(id_i-1)*K + id_j]])
      ) %>% pivot_longer(cols = 2:3) %>% dplyr::rename(models = "name", density = "value") %>% ggplot() +
        geom_area(aes(x = x, y = density, fill = models), alpha = 0.6, position = "identity") +
        ggtitle(paste0(id_i, "->", id_j, ", mixed")) + theme(legend.position="bottom")
    }
  }

  plist <- list(
    p_common, ps[[1]], pms[[1]],
    p_common, ps[[2]], pms[[2]],
    p_common, ps[[3]], pms[[3]],
    p_common, ps[[4]], pms[[4]]
  )
  ggsave(gridExtra::grid.arrange(grobs = plist, ncol=3, nrow = 4),
         file = file.path(fig_dir, paste0(filename, ".pdf")),
         width = 6, height = 8)
  return(res = list(RMISE = RMISE,
                    pws_cvg = pws_cvg,
                    cw = cw, is = is))
}
