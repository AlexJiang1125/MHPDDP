# do the same thing, but via sampling
calc_rmse_vi <- function(alpha0_vi, alphas_vi, cas_vi, cbs_vi, das_vi, dbs_vi,
                         ca0_vi, cb0_vi, da0_vi, db0_vi, eps_vi, type = "RANDOM",
                         fig_dir, ids, EXAMPLE_ID = EXAMPLE_ID) {
  library(EnvStats)
  res_list <- get_DPpars(EXAMPLE_ID = EXAMPLE_ID, samepriors = samepriors, epsilon_true = eps_true)
  list2env(res_list, envir = globalenv())
  res_list <- get_HPpars(HP_ID = HP_ID)
  list2env(res_list, envir = globalenv())
  # number of samples
  B_iter <- 5000
  # number of Beta mixture components
  L <- length(ca0_vi)
  # number of samples generated from the posterior variational distribution
  if (EXAMPLE_ID %in% c(6,8)) {
    PLOT_START <- 0.05
  } else {
    PLOT_START <- 0.01
  }
  x_plot <- seq(PLOT_START, T_kernel-PLOT_START, by = PLOT_START)
  n_xplot <- length(x_plot)
  # generate samples from the variational distribution
  # p0m <- alpha0_vi/sum(alpha0_vi)
  p0m <- rdirichlet(B_iter, alpha = alpha0_vi)
  a0m <- b0m <- matrix(0, nrow = B_iter, ncol = L)
  for (i in 1:B_iter) {
    a0m[i,] <- rgamma(L, ca0_vi, da0_vi)
    b0m[i,] <- rgamma(L, cb0_vi, db0_vi)
  }
  psm <- asm <- bsm <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      psm[[(id_i - 1)*K + id_j]] <- rdirichlet(B_iter, alpha = alphas_vi[[(id_i - 1)*K + id_j]]) #alphas_vi[[(id_i - 1)*K + id_j]]/sum(alphas_vi[[(id_i - 1)*K + id_j]])
      asm[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter, ncol = L)
      bsm[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter, ncol = L)
      for (i in 1:B_iter) {
        asm[[(id_i - 1)*K + id_j]][i,] <- rgamma(L, cas_vi[[(id_i - 1)*K + id_j]], das_vi[[(id_i - 1)*K + id_j]])
        bsm[[(id_i - 1)*K + id_j]][i,] <- rgamma(L, cbs_vi[[(id_i - 1)*K + id_j]], dbs_vi[[(id_i - 1)*K + id_j]])
      }
      #asm[[(id_i - 1)*K + id_j]] <- cas_vi[[(id_i - 1)*K + id_j]]/das_vi[[(id_i - 1)*K + id_j]]
      #bsm[[(id_i - 1)*K + id_j]] <- cbs_vi[[(id_i - 1)*K + id_j]]/dbs_vi[[(id_i - 1)*K + id_j]]
    }
  }
  if (type == "RANDOM") {
    epm <- rbeta(B_iter, eps_vi[1], eps_vi[2])  #eps_vi[1]/sum(eps_vi)
  } else if (type == "COMMON") {
    epm <- rep(1, B_iter)
  } else {
    epm <- rep(0, B_iter)
  }
  
  ## calculate true density on a grid
  # 1. common component part
  pos_p0s_true <- pos_f0s_true <- NA
  if (EXAMPLE_ID == 5) {
    for (j in 1:length(x_plot)) {
      pos_p0s_true[j] <- eps_tri0*dtri(x_plot[j], min = tri_a0[1], max = tri_c0[1], mode = tri_b0[1]) +
        (1-eps_tri0)*dtri(x_plot[j], min = tri_a0[2], max = tri_c0[2], mode = tri_b0[2])
    }
  } else if (EXAMPLE_ID %in% c(6,8)) {
    for (j in 1:length(x_plot)) {
      pos_p0s_true[j] <- dexp(x_plot[j], exp_0)
      pos_f0s_true[j] <- pexp(x_plot[j], exp_0)
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
      pos_p0s_true[j] <- dbeta(x_plot[j], a0s, b0s)
      pos_f0s_true[j] <- pbeta(x_plot[j], a0s, b0s)
    }
  }
  # 2. heterogenous components part
  pos_ps_true <- pos_fs_true <- list()
  pos_ps_true_mixed <- pos_fs_true_mixed <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos <- (id_i - 1)*K + id_j
      temp <- temp_f <- NA
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
          temp_f[j] <- pexp(x_plot[j], exp_s[[id_i]][[id_j]])
        } else {
          # for examples 4 (single mixture) and 7 (DP mixture)
          for (id_l in 1:length_pmixtures) {
            if (id_l == 1) {
              temp[j] <- ps[[id_i]][[id_j]][1]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][1], bs[[id_i]][[id_j]][1], 0, T_kernel)
              temp_f[j] <- ps[[id_i]][[id_j]][1]*pnsbeta(x_plot[j], as[[id_i]][[id_j]][1], bs[[id_i]][[id_j]][1], 0, T_kernel)
            } else {
              temp[j] <- temp[j] + ps[[id_i]][[id_j]][id_l]*dnsbeta(x_plot[j], as[[id_i]][[id_j]][id_l], bs[[id_i]][[id_j]][id_l], 0, T_kernel)
              temp_f[j] <- temp[j] + ps[[id_i]][[id_j]][id_l]*pnsbeta(x_plot[j], as[[id_i]][[id_j]][id_l], bs[[id_i]][[id_j]][id_l], 0, T_kernel)
            }
          }
        }
      }
      pos_ps_true[[(id_i - 1)*K + id_j]] <- temp
      pos_ps_true_mixed[[(id_i - 1)*K + id_j]] <- pos_p0s_true*(eps_true) + pos_ps_true[[(id_i - 1)*K + id_j]]*(1 - eps_true)
      pos_fs_true[[(id_i - 1)*K + id_j]] <- temp_f
      pos_fs_true_mixed[[(id_i - 1)*K + id_j]] <- pos_f0s_true*(eps_true) + pos_fs_true[[(id_i - 1)*K + id_j]]*(1 - eps_true)
    }
  }
  
  
  ## calculate posterior results
  pos_p0s <- pos_f0s <- matrix(0, nrow = B_iter, ncol = length(x_plot))
  for (j in 1:length(x_plot)) {
    pos_p0s[,j] <- eval_mixture(x_plot[j],p0m,a0m,b0m, T_kernel = T_kernel)
    pos_f0s[,j] <- eval_mixture(x_plot[j],p0m,a0m,b0m,type = "p", T_kernel = T_kernel)
  }
  
  # idio part
  pos_ps <- list()
  # mixed part
  pos_ps_mixed <- list()
  
  pos_fs <- pos_fs_mixed <- list()
  
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos_ps[[(id_i - 1)*K + id_j]] <- pos_fs[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter, ncol = length(x_plot))
      temp <- NA
      for (j in 1:length(x_plot)) {
        pos_ps[[(id_i - 1)*K + id_j]][,j] <- eval_mixture(x_plot[j], psm[[(id_i - 1)*K + id_j]],
                                                          asm[[(id_i - 1)*K + id_j]],
                                                          bsm[[(id_i - 1)*K + id_j]], T_kernel = T_kernel)
        pos_fs[[(id_i - 1)*K + id_j]][,j] <- eval_mixture(x_plot[j], psm[[(id_i - 1)*K + id_j]],
                                                          asm[[(id_i - 1)*K + id_j]],
                                                          bsm[[(id_i - 1)*K + id_j]], type = "p")
        
      }
      pos_ps_mixed[[(id_i - 1)*K + id_j]] <-
        pos_p0s * matrix(rep(epm, n_xplot), ncol = n_xplot) +
        pos_ps[[(id_i - 1)*K + id_j]] * matrix(rep((1-epm), n_xplot), ncol = n_xplot)
      pos_fs_mixed[[(id_i - 1)*K + id_j]] <-
        pos_f0s * matrix(rep(epm, n_xplot), ncol = n_xplot) +
        pos_fs[[(id_i - 1)*K + id_j]] * matrix(rep((1-epm), n_xplot), ncol = n_xplot)
      
    }
  }
  
  plot(pos_ps_true_mixed[[2]])
  plot(colMeans(pos_ps_mixed[[2]]))
  
  # calculate
  rmse_list <- list()
  wass_list <- list()
  aocs <- NA
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      rmse_list[[(id_i - 1)*K + id_j]] <- PLOT_START*(sweep(x = pos_ps_mixed[[(id_i - 1)*K + id_j]],
                                                            MARGIN = 2,
                                                            STATS = pos_ps_true_mixed[[(id_i - 1)*K + id_j]]))^2
      wass_list[[(id_i - 1)*K + id_j]] <- abs(sweep(x = pos_fs_mixed[[(id_i - 1)*K + id_j]],
                                                    MARGIN = 2,
                                                    STATS = pos_fs_true_mixed[[(id_i - 1)*K + id_j]]))
      aocs[(id_i - 1)*K + id_j] <- sum(colMeans(pos_ps_mixed[[(id_i - 1)*K + id_j]])*PLOT_START)
    }
  }
  rmse_agg <- lapply(rmse_list, FUN = "rowSums")
  wass_agg <- lapply(wass_list, FUN = "colMeans")
  RMISE <- NA
  WASS <- NA
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      RMISE[(id_i - 1)*K + id_j] <- mean(sqrt(rmse_agg[[(id_i - 1)*K + id_j]]))
      WASS[(id_i - 1)*K + id_j] <- sum(wass_agg[[(id_i - 1)*K + id_j]])*PLOT_START
    }
  }
  
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
  
  res = list(RMISE = RMISE,
             pws_cvg = pws_cvg,
             cw = cw, is = is, aoc = mean(aocs))
  
  return(res)
}

get_ELBO <- function() {
  # evaluate ELBO
  
  elbo_mu <- sum(-eta_mus1/eta_mus2*(T_all*split_ratio) + (digamma(eta_mus1) - log(eta_mus2))*(new_mus1 - 1))
  elbo_alpha <- sum(-eta_alphas1/eta_alphas2*(new_alphas2) + (new_alphas1 - 1)*(digamma(eta_alphas1) - log(eta_alphas2)))
  
  # p and eps
  
  elbo_a1 <- elbo_a2 <- elbo_b1 <- elbo_b2 <- elbo_B <- 0
  
  # update hyperpar for Dirichlet (ps)
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      id_ij <- (id_i - 1)*K + id_j
      
      ids_dim <- ids[[(id_i - 1)*K + id_j]]
      y_dim <- (ids_dim$time_end - ids_dim$time_start)
      B_dim <- B[cbind(ids_dim$Var2_new, ids_dim$Var1_new)]
      if (id_i == 1 && id_j == 1) {
        y0 <- y_dim
        B0 <- B_dim
        Z0_vi <- Zs_vi[[(id_i - 1)*K + id_j]][1:L, ]
      } else {
        Z0_vi <- cbind(Z0_vi, Zs_vi[[(id_i - 1)*K + id_j]][1:L, ])
        y0 <- c(y0, y_dim)
        B0 <- c(B0, B_dim)
      }
      tmp <- Zs_vi[[(id_i - 1)*K + id_j]]
      if (id_i == 1 && id_j == 1) {
        tmp2 <- tmp
      } else {
        tmp2 <- cbind(tmp2,tmp)
      }
      elbo_a1 <- elbo_a1 + sum(apply(tmp[1:L + L,] * matrix(rep(B_dim, L), nrow = L, byrow = TRUE),
                                     MARGIN = 1,FUN = "sum") * us_bar[[id_ij]] * (
                                       digamma(us_bar[[id_ij]] + vs_bar[[id_ij]]) - digamma(us_bar[[id_ij]]) + vs_bar[[id_ij]] *
                                         psigamma(us_bar[[id_ij]] + vs_bar[[id_ij]], deriv = 1) * (qlbs[[id_ij]] - log(vs_bar[[id_ij]]))
                                     )*(digamma(cas_vi[[id_ij]]) - log(das_vi[[id_ij]])))
      elbo_a2 <- elbo_a2 - sum((cas_vi[[id_ij]]/das_vi[[id_ij]] - 1)*apply((Zs_vi[[id_ij]][(L+1):(2*L), ])*matrix(rep(log(y_dim/T_kernel)*B_dim, L), nrow = L, byrow = TRUE), MARGIN = 1, FUN = "sum"))
      
      elbo_b1 <- elbo_b1 + sum(apply(tmp[1:L + L,] * matrix(rep(B_dim, L), nrow = L, byrow = TRUE),
                                     MARGIN = 1,FUN = "sum") * vs_bar[[id_ij]] * (
                                       digamma(vs_bar[[id_ij]] + us_bar[[id_ij]]) - digamma(vs_bar[[id_ij]]) + us_bar[[id_ij]] *
                                         psigamma(vs_bar[[id_ij]] + us_bar[[id_ij]], deriv = 1) * (qlas[[id_ij]] - log(us_bar[[id_ij]]))
                                     )*(digamma(cbs_vi[[id_ij]]) - log(dbs_vi[[id_ij]])))
      elbo_b2 <- elbo_b2 - sum((cbs_vi[[id_ij]]/dbs_vi[[id_ij]] - 1)*apply((Zs_vi[[(id_i - 1)*K + id_j]][(L+1):(2*L), ])*matrix(rep(log(1 - y_dim/T_kernel)*B_dim, L), nrow = L, byrow = TRUE), MARGIN = 1, FUN = "sum"))
    }
  }
  elbo_a1 <- elbo_a1 + sum((digamma(ca0_vi) - log(da0_vi))*apply(Z0_vi*matrix(rep(B0,L), nrow = L, byrow = TRUE),
                                                                 MARGIN = 1, FUN = "sum")/split_ratio*u0_bar*(digamma(u0_bar + v0_bar) - digamma(u0_bar) + v0_bar*psigamma(u0_bar + v0_bar, deriv = 1)*(qlb0 - log(v0_bar))))
  elbo_a2 <- elbo_a2-sum((ca0_vi/da0_vi - 1)*apply(Z0_vi*matrix(rep((log(y0) - log(T_kernel))*B0, L), nrow = L, byrow = TRUE), MARGIN = 1, FUN = "sum"))
  
  elbo_b1 <- elbo_b1 + sum((digamma(cb0_vi) - log(db0_vi))*apply(Z0_vi*matrix(rep(B0,L), nrow = L, byrow = TRUE),
                                                                 MARGIN = 1, FUN = "sum")/split_ratio*u0_bar*(digamma(u0_bar + v0_bar) - digamma(u0_bar) + v0_bar*psigamma(u0_bar + v0_bar, deriv = 1)*(qla0 - log(v0_bar))))
  elbo_b2 <- elbo_b2-sum((cb0_vi/db0_vi - 1)*apply(Z0_vi*matrix(rep((log(T_kernel-y0) - log(T_kernel))*B0, L), nrow = L, byrow = TRUE), MARGIN = 1, FUN = "sum"))
  
  eps_part1 <- digamma(eps_vi[1]) - digamma(eps_vi[1] + eps_vi[2])
  eps_part2 <- digamma(eps_vi[2]) - digamma(eps_vi[1] + eps_vi[2])
  # update eps
  elbo_eps <- (sum(tmp2[1:L, ]*matrix(rep(B0,L), nrow = L, byrow = TRUE)) - 1)*(eps_part1) +
    (sum(tmp2[1:L+L, ]*matrix(rep(B0,L), nrow = L, byrow = TRUE)) - 1)*(eps_part2)
  # update p
  elbo_p <- sum((digamma(alpha0_vi) - sum(digamma(alpha0_vi)))*apply(tmp2[1:L, ]*matrix(rep(B0,L), nrow = L, byrow = TRUE), MARGIN = 1, FUN = "sum"))
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      elbo_p <- elbo_p + sum((digamma(alphas_vi[[(id_i - 1)*K + id_j]]) - sum(digamma(alphas_vi[[(id_i - 1)*K + id_j]])))*apply(tmp[(L+1):(2*L), ]*matrix(rep(B_dim,L), nrow = L, byrow = TRUE), MARGIN = 1, FUN = "sum"))
    }
  }
  # update B
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      id_ij <- (id_i - 1)*K + id_j
      ids_dim <- ids[[(id_i - 1)*K + id_j]]
      y_dim <- (ids_dim$time_end - ids_dim$time_start)
      B_dim <- B[cbind(ids_dim$Var2_new, ids_dim$Var1_new)]
      elbo_B <- elbo_B + sum(B_dim*(log(T_kernel) - log(y_dim) - log(T_kernel - y_dim)))
    }
  }
  
  elbo <- elbo_p + elbo_eps + elbo_b2 + elbo_b1 + elbo_a1 + elbo_a2 + elbo_B
  return(elbo)
}

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
