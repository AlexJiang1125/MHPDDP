calc_elbo <- function(elbo_parts = 5) {
  # split into five parts
  elbos <- NA
  for (id_elbo in 1:elbo_parts) {
    # generate random start and end dates
    # T_all <- length(y[[1]])
    time_start <- T_all * split_ratio *(id_elbo-1)
    time_end <- T_all * split_ratio * id_elbo
    # y is the overall
    out <- out_all[out_all$timestamp < time_end & out_all$timestamp >= time_start ,]
    out$timestamp <- out$timestamp - time_start
    id1s <- NA
    ids <- list()
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        tmp <- ids_all[[(id_i - 1) * K + id_j]]
        ids[[(id_i - 1) * K + id_j]] <- tmp[tmp$time_start >= time_start & tmp$time_end < time_end ,]
        ids[[(id_i - 1) * K + id_j]]$time_start <- ids[[(id_i - 1) * K + id_j]]$time_start - time_start
        ids[[(id_i - 1) * K + id_j]]$time_end <- ids[[(id_i - 1) * K + id_j]]$time_end - time_start
        id1s[(id_i - 1) * K + id_j] <- min(ids[[(id_i - 1) * K + id_j]]$Var1)
      }
    }
    # update id's
    id_min <- as.integer(rownames(out)[1])
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        ids[[(id_i - 1) * K + id_j]]$Var1_new <- ids[[(id_i - 1) * K + id_j]]$Var1 - id_min + 1
        ids[[(id_i - 1) * K + id_j]]$Var2_new <- ids[[(id_i - 1) * K + id_j]]$Var2 - id_min + 1

      }
    }

    out <- res$out
    ids <- res$ids
    N <- dim(out)[1]
    Ns <- sapply(ids, FUN = "dim")[1,]
    Ns_alpha <- as.vector(table(out$dim))

    ## update B
    ## update assignment
    # intermediate updates
    u0_bar <- ca0_vi/da0_vi
    v0_bar <- cb0_vi/db0_vi
    # E[ln(pi_i)]
    qlpi0 <- digamma(alpha0_vi) - digamma(sum(alpha0_vi))
    # E[ln(a)]
    qla0 <- digamma(ca0_vi) - log(da0_vi)
    # E[ln(b)]
    qlb0 <- digamma(cb0_vi) - log(db0_vi)
    # E[(ln(u) - ln(ubar))^2]
    qlusq0 <- (digamma(ca0_vi) - log(ca0_vi))^2 + psigamma(ca0_vi, deriv = 1)
    # E[(ln(v) - ln(vbar))^2]
    qlvsq0 <- (digamma(cb0_vi) - log(cb0_vi))^2 + psigamma(cb0_vi, deriv = 1)
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        us_bar[[(id_i - 1)*K + id_j]] <- cas_vi[[(id_i - 1)*K + id_j]]/das_vi[[(id_i - 1)*K + id_j]]
        vs_bar[[(id_i - 1)*K + id_j]] <- cbs_vi[[(id_i - 1)*K + id_j]]/dbs_vi[[(id_i - 1)*K + id_j]]
        qlpis[[(id_i - 1)*K + id_j]] <- digamma(alphas_vi[[(id_i - 1)*K + id_j]]) - digamma(sum(alphas_vi[[(id_i - 1)*K + id_j]]))
        qlas[[(id_i - 1)*K + id_j]] <- digamma(cas_vi[[(id_i - 1)*K + id_j]]) - log(das_vi[[(id_i - 1)*K + id_j]])
        qlbs[[(id_i - 1)*K + id_j]] <- digamma(cbs_vi[[(id_i - 1)*K + id_j]]) - log(dbs_vi[[(id_i - 1)*K + id_j]])
        qlusqs[[(id_i - 1)*K + id_j]] <- (digamma(cas_vi[[(id_i - 1)*K + id_j]]) - log(cas_vi[[(id_i - 1)*K + id_j]]))^2 +
          psigamma(cas_vi[[(id_i - 1)*K + id_j]], deriv = 1)
        qlvsqs[[(id_i - 1)*K + id_j]] <- (digamma(cbs_vi[[(id_i - 1)*K + id_j]]) - log(cbs_vi[[(id_i - 1)*K + id_j]]))^2 +
          psigamma(cbs_vi[[(id_i - 1)*K + id_j]], deriv = 1)
      }
    }
    Zs_vi <- list()
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        Zs_vi[[(id_i - 1)*K + id_j]] <-  t(rdirichlet(Ns[(id_i - 1)*K + id_j], alpha = rep(alpha_DP/L, 2*L))) #matrix(1/L, nrow = L, ncol = N)
      }
    }
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        ids_dim <- ids[[(id_i - 1)*K + id_j]]
        y_dim <- ids_dim$time_end - ids_dim$time_start
        # update common components
        for (l in 1:L) {
          if (borrow_mode == "RANDOM") {
            eps_part <- digamma(eps_vi[1]) - digamma(eps_vi[1] + eps_vi[2])
          } else if (borrow_mode == "COMMON") {
            eps_part <- 0
          } else {
            eps_part <- -Inf
          }
          Zs_vi[[(id_i - 1)*K + id_j]][l, ] <- qlpi0[l] + (u0_bar[l] - 1)*log(y_dim) + (v0_bar[l] - 1)*log(T_kernel - y_dim) -
            (u0_bar[l] + v0_bar[l] - 1)*log(T_kernel) -
            lbeta(u0_bar[l], v0_bar[l]) +
            u0_bar[l]*(digamma(u0_bar[l] + v0_bar[l]) - digamma(u0_bar[l]))*(qla0[l] - log(u0_bar[l])) +
            v0_bar[l]*(digamma(u0_bar[l] + v0_bar[l]) - digamma(v0_bar[l]))*(qlb0[l] - log(v0_bar[l])) +
            0.5*u0_bar[l]^2*(psigamma(u0_bar[l] + v0_bar[l], deriv = 1) - psigamma(u0_bar[l], deriv = 1))*qlusq0[l] +
            0.5*v0_bar[l]^2*(psigamma(u0_bar[l] + v0_bar[l], deriv = 1) - psigamma(v0_bar[l], deriv = 1))*qlvsq0[l] +
            u0_bar[l]*v0_bar[l]*psigamma(u0_bar[l] + v0_bar[l], deriv = 1)*(qla0[l] - log(u0_bar[l]))*(qlb0[l] - log(v0_bar[l])) +
            eps_part
          if (borrow_mode == "RANDOM") {
            eps_part <- digamma(eps_vi[2]) - digamma(eps_vi[1] + eps_vi[2])
          } else if (borrow_mode == "COMMON") {
            eps_part <- -Inf
          } else {
            eps_part <- 0
          }

          Zs_vi[[(id_i - 1)*K + id_j]][l+L, ] <- qlpis[[(id_i - 1)*K + id_j]][l] + (us_bar[[(id_i - 1)*K + id_j]][l] - 1)*log(y_dim) +
            (vs_bar[[(id_i - 1)*K + id_j]][l] - 1)*log(T_kernel - y_dim) - (us_bar[[(id_i - 1)*K + id_j]][l] + vs_bar[[(id_i - 1)*K + id_j]][l] - 1)*log(T_kernel) -
            lbeta(us_bar[[(id_i - 1)*K + id_j]][l], vs_bar[[(id_i - 1)*K + id_j]][l]) +
            us_bar[[(id_i - 1)*K + id_j]][l]*(digamma(us_bar[[(id_i - 1)*K + id_j]][l] + vs_bar[[(id_i - 1)*K + id_j]][l]) - digamma(us_bar[[(id_i - 1)*K + id_j]][l]))*(qlas[[(id_i - 1)*K + id_j]][l] - log(us_bar[[(id_i - 1)*K + id_j]][l])) +
            vs_bar[[(id_i - 1)*K + id_j]][l]*(digamma(us_bar[[(id_i - 1)*K + id_j]][l] + vs_bar[[(id_i - 1)*K + id_j]][l]) - digamma(vs_bar[[(id_i - 1)*K + id_j]][l]))*(qlbs[[(id_i - 1)*K + id_j]][l] - log(vs_bar[[(id_i - 1)*K + id_j]][l])) +
            0.5*us_bar[[(id_i - 1)*K + id_j]][l]^2*(psigamma(us_bar[[(id_i - 1)*K + id_j]][l] + vs_bar[[(id_i - 1)*K + id_j]][l], deriv = 1) - psigamma(us_bar[[(id_i - 1)*K + id_j]][l], deriv = 1))*qlusqs[[(id_i - 1)*K + id_j]][l] +
            0.5*vs_bar[[(id_i - 1)*K + id_j]][l]^2*(psigamma(us_bar[[(id_i - 1)*K + id_j]][l] + vs_bar[[(id_i - 1)*K + id_j]][l], deriv = 1) - psigamma(vs_bar[[(id_i - 1)*K + id_j]][l], deriv = 1))*qlvsqs[[(id_i - 1)*K + id_j]][l] +
            us_bar[[(id_i - 1)*K + id_j]][l]*vs_bar[[(id_i - 1)*K + id_j]][l]*psigamma(us_bar[[(id_i - 1)*K + id_j]][l] + vs_bar[[(id_i - 1)*K + id_j]][l], deriv = 1)*(qlas[[(id_i - 1)*K + id_j]][l] - log(us_bar[[(id_i - 1)*K + id_j]][l]))*(qlbs[[(id_i - 1)*K + id_j]][l] - log(vs_bar[[(id_i - 1)*K + id_j]][l])) +
            eps_part
        }
      }
    }

    # convert to r_nl
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        Zs_vi[[(id_i - 1)*K + id_j]] <- apply(Zs_vi[[(id_i - 1)*K + id_j]], MARGIN = 2, FUN = "exp_norm")
      }
    }

    cas_vi_mat <- matrix(unlist(cas_vi), nrow = K^2, byrow = TRUE)
    cbs_vi_mat <- matrix(unlist(cbs_vi), nrow = K^2, byrow = TRUE)
    das_vi_mat <- matrix(unlist(das_vi), nrow = K^2, byrow = TRUE)
    dbs_vi_mat <- matrix(unlist(dbs_vi), nrow = K^2, byrow = TRUE)

    # a vector of E[lbeta(a0, b0)], with L elements
    lbeta0s <- exp_lbeta(ca0_vi, cb0_vi, da0_vi, db0_vi)
    a0_bar1 <- ca0_vi/da0_vi - 1
    b0_bar1 <- cb0_vi/db0_vi - 1
    if (split_ratio == 1) {
      B <- Matrix::Matrix(nrow = N, ncol = N, data = 0, sparse = TRUE,
                          doDiag = FALSE)
      B[1,1] <- 1
      #se <- as(Bsparse, "dgTMatrix")
    } else {
      B <- matrix(0, N, N)
      B[1, 1] <- 1
    }
    diag(B) <- digamma(eta_mus1[out$dim]) - log(eta_mus2[out$dim])
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        lbetas <- exp_lbeta(cas_vi_mat[(id_i - 1)*K + id_j,],
                            cbs_vi_mat[(id_i - 1)*K + id_j,],
                            das_vi_mat[(id_i - 1)*K + id_j,],
                            dbs_vi_mat[(id_i - 1)*K + id_j,])
        a_bar1 <- cas_vi_mat[(id_i - 1)*K + id_j,]/das_vi_mat[(id_i - 1)*K + id_j,] - 1
        a_bar1_m <- matrix(rep(a_bar1, Ns[(id_i - 1)*K + id_j]), nrow = L, byrow = FALSE)
        a0_bar1_m <- matrix(rep(a0_bar1, Ns[(id_i - 1)*K + id_j]), nrow = L, byrow = FALSE)
        b_bar1 <- cbs_vi_mat[(id_i - 1)*K + id_j,]/dbs_vi_mat[(id_i - 1)*K + id_j,] - 1
        b_bar1_m <- matrix(rep(b_bar1, Ns[(id_i - 1)*K + id_j]), nrow = L, byrow = FALSE)
        b0_bar1_m <- matrix(rep(b0_bar1, Ns[(id_i - 1)*K + id_j]), nrow = L, byrow = FALSE)
        ids_dim <- ids[[(id_i - 1)*K + id_j]]
        y_dim <- ids_dim$time_end - ids_dim$time_start
        y_dim_m <- matrix(rep(log(y_dim) - log(T_kernel), each = L), nrow = L, byrow = FALSE)
        yneg_dim_m <- matrix(rep(log(T_kernel-y_dim) - log(T_kernel), each = L), nrow = L, byrow = FALSE)
        T_kernel_dim_m <- matrix(rep(log(rep(T_kernel, length(y_dim))), each = L), nrow = L, byrow = FALSE)
        lbetas_m <- matrix(rep(lbetas, Ns[(id_i - 1)*K + id_j]), nrow = L, byrow = FALSE)
        lbeta0s_m <- matrix(rep(lbeta0s, Ns[(id_i - 1)*K + id_j]), nrow = L, byrow = FALSE)
        part1 <- digamma(eta_alphas1[(id_i - 1)*K + id_j]) - log(eta_alphas2[(id_i - 1)*K + id_j]) - log(T_kernel) # alpha parts
        part2 <- apply(Zs_vi[[(id_i - 1)*K + id_j]][1:L+L,]*(a_bar1_m*y_dim_m + b_bar1_m*yneg_dim_m + lbetas_m), FUN = "sum", 2)
        part3 <- apply(Zs_vi[[(id_i - 1)*K + id_j]][1:L,]*(a0_bar1_m*y_dim_m + b0_bar1_m*yneg_dim_m + lbeta0s_m), FUN = "sum", 2)
        B[cbind(ids_dim$Var2_new, ids_dim$Var1_new)] <- B[cbind(ids_dim$Var2_new, ids_dim$Var1_new)] + part1 + part2 + part3

      }
    }

    # B <- do.call("rbind", rowapply_simple_triplet_matrix(as.simple_triplet_matrix(B),
    #                                                      FUN = "exp_norm_rmzero"))

    if (split_ratio == 1) {
      Blist<- lapply(listRows(B), FUN = "exp_norm_rmzero")
      nlen <- unlist(lapply(Blist, FUN = "length"))
      sparse_idi <- rep(1:dim(B)[1], times = nlen)
      id_end <- 1:dim(B)[1]
      id_start <- id_end - nlen + 1
      id_mt <- cbind(id_start, id_end)
      id_list <- lapply(seq_len(nrow(id_mt)), function(i) id_mt[i,])
      sparse_idj <- unlist(lapply(id_list, FUN = "my_seqint"))
      B <- sparseMatrix(i = sparse_idi,
                        j = sparse_idj,
                        x= unlist(Blist))
      rm(Blist)
    } else {
      B <- t(apply(B, FUN = "exp_norm_rmzero", MARGIN = 1))
    }




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
    elbos[id_elbo] <- elbo_p + elbo_eps + elbo_b2 + elbo_b1 + elbo_a1 + elbo_a2 + elbo_B

  }
  return(sum(elbos))
}
