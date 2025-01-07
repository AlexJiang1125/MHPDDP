MHPDDP_mcmc <- function(data = df,
                        K = 2, e_alpha = 2, f_alpha = 4,
                        T_all = 15000, a_mu = 1, b_mu = 100,
                        B_iter = 10001, delta_mh = 0.015,
                        T_kernel = 1, c_a = 2, c_b = 3.5,
                        c_a0 = 2, c_b0 = 2, d_a0 = 2, d_b0 = 2,
                        seed = 123, borrow_mode = c("COMMON", "IDIO", "RANDOM")
                        ) {

  # read in the data
  y <- df$timestamp
  N <- length(y)

  # DP related hyperparameters
  # true parameters
  # dimensions for the MHP process

  start_time <- Sys.time()
  # HP parameters

  # common-idiosyncratic mixture weight
  epsilon <- rbeta(1, 1, 1)

  # a change in data structure
  Bvec <- rep(0, N)
  # assign values to true latent structure
  Bvec[1] <- 1
  for (i in 2:N) {
    if (out$parentdim[i] == 0) {
      Bvec[i] <- i
    } else {
      parentid <- out$parent[i]
      parentrowid <- which(out$id == parentid)
      Bvec[i] <- parentrowid
    }
  }

  # a compact way to store the branching structure matrix
  # only look at the closest 200 observations
  col_max <- 200
  B <- B_true <- matrix(0, N, col_max)
  #temp <- 0
  ydelta <- NA
  for (i in 1:N) {
    if (Bvec[i] == i) {
      B_true[i,1] <- 1
      ydelta[i] <- 0
    } else {
      #temp <- max(temp, i - Bvec[i])
      x <- min(i - Bvec[i] + 1, col_max)
      B_true[i,x] <- 1
      ydelta[i] <- y[i] - y[Bvec[i]]
    }
  }

  for (i in 1:N) {
    if (Bvec[i] == i) {
      B[i,1] <- 1
    } else {
      x <- min(i - Bvec[i] + 1, col_max)
      if (ydelta[i] >= T_kernel) {
        B[i,1] <- 1
      } else {
        B[i,x] <- 1
      }
    }
  }

  # get transition pairs
  ids <- list()
  pos <- 1
  for (j in 1:K) {
    for (k in 1:K) {
      id_j <- which(out$dim == j)
      id_k <- which(out$dim == k)
      id_jk <- expand.grid(id_j, id_k)
      id_jk <- id_jk[id_jk$Var2 - id_jk$Var1 > 0, ]
      id_jk <- id_jk[id_jk$Var2 - id_jk$Var1 < col_max, ]
      #id_jk <- id_jk[out$parentdim[id_jk$Var1] != 0 , ]
      ids[[pos]] <- id_jk
      pos <- pos + 1
    }
  }

  out_tail <- out[out$timestamp < T_all - 0, ]
  ys_truncated <- list()
  for (id_k in 1:K) {
    ys_truncated[[id_k]] <- out_tail[out_tail$dim == id_k,]
  }

  # HP MCMC parameters
  alphas_mcmc  <- list()#output$alpha_sgld[c(1:count)]
  mus_mcmc  <- list()#output$beta_sgld[c(1:count)]

  a0_mcmc <- matrix(0, nrow = B_iter, ncol = L)
  b0_mcmc <- matrix(0, nrow = B_iter, ncol = L)
  p0_mcmc <- matrix(0, nrow = B_iter, ncol = L)
  as_mcmc <- list()
  bs_mcmc <- list()
  p_mcmc <- list()
  Z_mcmc <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      as_mcmc[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter, ncol = L)
      bs_mcmc[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter, ncol = L)
      p_mcmc[[(id_i - 1)*K + id_j]] <- matrix(0, nrow = B_iter, ncol = L)
      Z_mcmc[[(id_i - 1)*K + id_j]] <- rep(0, N)
    }
  }

  set.seed(seed) # for 3 mixture model, covers truth
  alphas_mcmc  <- list()#output$alpha_sgld[c(1:count)]
  mus_mcmc  <- list()#output$beta_sgld[c(1:count)]

  alpha <- matrix(rgamma(K*K, e_alpha, f_alpha), nrow = K)
  mu <- rgamma(K, a_mu, b_mu)

  B_immi <- 1
  B_overall <- 1

  a0_mcmc[1,] <- rgamma(L, c_a, d_a)
  b0_mcmc[1,] <- rgamma(L, c_b, d_b)
  p0_mcmc[1,] <- rdirichlet(1, alpha = rep(alpha_DP/L, L))[1,]

  for (id_i in 1:K) {
    for (id_j in 1:K) {
      as_mcmc[[(id_i - 1)*K + id_j]][1,] <- rgamma(L, c_a, d_a)
      bs_mcmc[[(id_i - 1)*K + id_j]][1,] <- rgamma(L, c_b, d_b)
      p_mcmc[[(id_i - 1)*K + id_j]][1,] <- rdirichlet(1, alpha = rep(alpha_DP/L, L))[1,]
    }
  }

  as_array <- simplify2array(as_mcmc)
  bs_array <- simplify2array(bs_mcmc)
  ps_array <- simplify2array(p_mcmc)

  for (ii in 2:N) {
    ptemp <- rep(0,min(col_max, ii))
    ptemp[1] <- log(mu[out$dim[ii]])
    id_start <- max(1, ii - col_max + 1)
    id_end <- ii - 1
    # find how many of them will be in the nonzero probability region
    # that's excluding y[ii] itself
    ids_nonzero <- sum((y[ii] - y[id_end:id_start]) < T_kernel)
    if (ids_nonzero == 0) {
      B[ii,1] <- 1
    } else {
      id_start_new <- ii - ids_nonzero
      y_evals <- y[ii] - y[id_end:id_start_new]
      ptemp[2:(ids_nonzero + 1)] <- log(alpha[out$dim[id_end:id_start_new], out$dim[ii]])
      # OMGGGG
      ids_pos <- (out$dim[id_end:id_start_new] - 1) * K + out$dim[ii]
      vtemp <- rep(0, length(y_evals))
      for (id_m in 1:L) {
        vtemp <- vtemp + dnsbeta(y_evals, as_array[1, id_m, ids_pos],
                                 bs_array[1, id_m, ids_pos],0, T_kernel)*ps_array[1, id_m, ids_pos]*(1-epsilon_mcmc[1]) +
          dnsbeta(y_evals, a0_mcmc[1, id_m] , b0_mcmc[1, id_m],0, T_kernel)*p0_mcmc[1, id_m]*(epsilon_mcmc[1])
      }
      ptemp[2:(ids_nonzero + 1)] <- ptemp[2:(ids_nonzero + 1)] + log(vtemp)
      ptemp <- ptemp[1:(ids_nonzero + 1)]
      res_temp <- exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp)))
      sampled <- sample(1:length(res_temp), size = 1, prob = res_temp)
      B[ii,] <- 0
      B[ii,sampled] <- 1#sample(1:min(col_max,ii), size = 1, prob = res_temp)
    }
  }

  for (id_i in 1:K) {
    for (id_j in 1:K) {
      pos <- (id_i - 1)*K + id_j
      ntrans <- sum(B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)])
      if (STORE_Z) {
        Z_mcmc[[(id_i - 1)*K + id_j]][1,] <- c(sample(1:L, size = ntrans, prob = p_mcmc[[(id_i - 1)*K + id_j]][1,], replace = TRUE), rep(0, N-ntrans))
      } else {
        Z_mcmc[[(id_i - 1)*K + id_j]] <- c(sample(1:(2*L), size = ntrans, prob = c(p0_mcmc[1,]*epsilon, p_mcmc[[(id_i - 1)*K + id_j]][1,]*(1-epsilon)), replace = TRUE), rep(0, N-ntrans))
      }

    }
  }

  for (i in 1:(B_iter-1)) {
    # update HP parts
    ntrans_all <- NA
    for (id_i in 1:K) {
      ndim <- sum(out$dim == id_i)
      for (id_j in 1:K) {
        pos <- (id_i - 1) * K + id_j
        ntrans <- sum(B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)])
        ntrans_all[pos] <- ntrans
        dists_B <- B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)]
        dists_y <- y[ids[[pos]]$Var2] - y[ids[[pos]]$Var1]
        dists <- sum(dists_B * dists_y)
        alpha_compensator <- 0
        for (j in 1:L) {
          alpha_compensator <- alpha_compensator + sum(pnsbeta(q = T_all - ys_truncated[[id_i]]$timestamp,
                                                               shape1 = as_mcmc[[(id_i-1)*K + id_j]][i,j],
                                                               shape2 = bs_mcmc[[(id_i-1)*K + id_j]][i,j],
                                                               min = 0, max = T_kernel)*p_mcmc[[(id_i-1)*K + id_j]][i,j])
        }
        truncation_correction <- alpha[id_i, id_j]*alpha_compensator
        alpha[id_i, id_j] <- rgamma(1, e_alpha + ntrans, f_alpha + alpha_compensator)
      }
      mu[id_i] <- rgamma(1, shape = sum(B[out$dim == id_i,1]) + a_mu, rate = T_all + b_mu)
    }
    # update Z
    # update common parts
    Z_all_mcmc <- c(Z_mcmc[[1]], Z_mcmc[[2]], Z_mcmc[[3]], Z_mcmc[[4]])

    counts0 <- as.vector(table(factor(Z_all_mcmc, levels = 1:L)))
    p0_mcmc[i+1, ] <- rdirichlet(1, alpha = alpha_DP/L + counts0)[1,]
    for (j in 1:L) {
      # joint proposal of a,b (on the log scale)
      a_star <- exp(rnorm(1, log(a0_mcmc[i,j]), sd = delta_mh))
      b_star <- exp(rnorm(1, log(b0_mcmc[i,j]), sd = delta_mh))
      # calculate un-normalized posterior of a,b
      # get all interarrival times between id_i and id_j
      # also those interarrivals are assigned to kernel j in the common comp.
      n_mix <- counts0[j]
      y_mix <- NA
      for (id_i in 1:K) {
        for (id_j in 1:K) {
          # get the interarrival times
          pos <- (id_i - 1) * K + id_j
          dists_B <- B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)]
          dists_y <- y[ids[[pos]]$Var2] - y[ids[[pos]]$Var1]
          interarr <- dists_B * dists_y
          interarr <- interarr[interarr > 0]
          tp <- Z_mcmc[[pos]]
          tp <- tp[tp > 0]
          if (length(interarr[tp == j]) > 0) {
            y_mix <- c(y_mix, interarr[tp == j])
          }
        }
      }
      y_mix <- y_mix[-1]
      ratio_part1 <- n_mix*(lbeta(a0_mcmc[i,j], b0_mcmc[i,j]) - lbeta(a_star, b_star) + log(T_kernel)*(a0_mcmc[i,j] + b0_mcmc[i,j] - a_star - b_star))
      ratio_part2 <- (a_star - a0_mcmc[i,j])*sum(log(y_mix)) + (b_star - b0_mcmc[i,j])*sum(log(T_kernel - y_mix))
      ratio_part3 <- (c_a0)*(log(a_star) - log(a0_mcmc[i,j])) + (c_b0)*(log(b_star) - log(b0_mcmc[i,j]))
      ratio_part4 <- -d_a0*(a_star - a0_mcmc[i,j]) - d_b0*(b_star - b0_mcmc[i,j])
      ratio <- ratio_part1 + ratio_part2 + ratio_part3 + ratio_part4
      if (runif(1, 0, 1) < exp(ratio)) {
        a0_mcmc[i+1, j] <- a_star
        b0_mcmc[i+1, j] <- b_star
      } else {
        a0_mcmc[i+1, j] <- a0_mcmc[i, j]
        b0_mcmc[i+1, j] <- b0_mcmc[i, j]
      }
    }
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        pos <- (id_i - 1) * K + id_j
        dists_B <- B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)]
        dists_y <- y[ids[[pos]]$Var2] - y[ids[[pos]]$Var1]
        interarr <- dists_B * dists_y
        y_DP <- interarr[interarr > 0]
        n_DP <- length(y_DP)
        counts <- as.vector(table(factor(Z_mcmc[[(id_i-1)*K + id_j]], levels = 1:L+L)))
        # update B
        p_mcmc[[(id_i-1)*K + id_j]][i+1, ] <- rdirichlet(1, alpha = alpha_DP/L + counts)[1,]
        # update atoms
        for (j in 1:L) {
          a_star <- exp(rnorm(1, log(as_mcmc[[(id_i-1)*K + id_j]][i,j]), sd = delta_mh))
          b_star <- exp(rnorm(1, log(bs_mcmc[[(id_i-1)*K + id_j]][i,j]), sd = delta_mh))
          # calculate un-normalized posterior of a,b
          counts <- as.vector(table(factor(Z_mcmc[[(id_i-1)*K + id_j]], levels = 1:L+L)))
          n_mix <- counts[j]
          ytemp <- y_DP
          y_mix <- ytemp[Z_mcmc[[(id_i-1)*K + id_j]] == j+L]
          ratio_part1 <- n_mix*(lbeta(as_mcmc[[(id_i-1)*K + id_j]][i,j], bs_mcmc[[(id_i-1)*K + id_j]][i,j]) - lbeta(a_star, b_star) + log(T_kernel)*(as_mcmc[[(id_i-1)*K + id_j]][i,j] + bs_mcmc[[(id_i-1)*K + id_j]][i,j] - a_star - b_star))
          ratio_part2 <- (a_star - as_mcmc[[(id_i-1)*K + id_j]][i,j])*sum(log(y_mix)) + (b_star - bs_mcmc[[(id_i-1)*K + id_j]][i,j])*sum(log(T_kernel - y_mix))
          ratio_part3 <- (c_a)*(log(a_star) - log(as_mcmc[[(id_i-1)*K + id_j]][i,j])) + (c_b)*(log(b_star) - log(bs_mcmc[[(id_i-1)*K + id_j]][i,j]))
          ratio_part4 <- -d_a*(a_star - as_mcmc[[(id_i-1)*K + id_j]][i,j]) - d_b*(b_star - bs_mcmc[[(id_i-1)*K + id_j]][i,j])
          ratio <- ratio_part1 + ratio_part2 + ratio_part3 + ratio_part4
          if (runif(1, 0, 1) < exp(ratio)) {
            as_mcmc[[(id_i-1)*K + id_j]][i+1, j] <- a_star
            bs_mcmc[[(id_i-1)*K + id_j]][i+1, j] <- b_star
          } else {
            as_mcmc[[(id_i-1)*K + id_j]][i+1, j] <- as_mcmc[[(id_i-1)*K + id_j]][i, j]
            bs_mcmc[[(id_i-1)*K + id_j]][i+1, j] <- bs_mcmc[[(id_i-1)*K + id_j]][i, j]
          }
        }
      }
    }


    # update epsilon
    if (borrow_mode == "COMMON") {
      epsilon_mcmc[i+1] <- 1
    } else if (borrow_mode == "IDIO") {
      epsilon_mcmc[i+1] <- 0
    } else {
      # count observations assigned to common components
      n_common <- n_idio <- 0
      for (id_i in 1:K) {
        for (id_j in 1:K) {
          n_common <- n_common + sum(Z_mcmc[[(id_i - 1)*K + id_j]] <= L) - sum(Z_mcmc[[(id_i - 1)*K + id_j]] == 0)
          n_idio <- n_idio + sum(Z_mcmc[[(id_i - 1)*K + id_j]] > L)
        }
      }
      epsilon_mcmc[i+1] <- rbeta(1, n_common + 1, n_idio + 1)
    }
    # update B
    as_array <- simplify2array(as_mcmc)
    bs_array <- simplify2array(bs_mcmc)
    ps_array <- simplify2array(p_mcmc)
    for (ii in 2:N) {
      ptemp <- rep(0,min(col_max, ii))
      ptemp[1] <- log(mu[out$dim[ii]])
      id_start <- max(1, ii - col_max + 1)
      id_end <- ii - 1
      # find how many of them will be in the nonzero probability region
      # that's excluding y[ii] itself
      ids_nonzero <- sum((y[ii] - y[id_end:id_start]) < T_kernel)

      # if all of them are outside, assign it to itself
      if (ids_nonzero == 0) {
        B[ii,1] <- 1
      } else {
        id_start_new <- ii - ids_nonzero
        y_evals <- y[ii] - y[id_end:id_start_new]
        ptemp[2:(ids_nonzero + 1)] <- log(alpha[out$dim[id_end:id_start_new], out$dim[ii]])
        # OMGGGG
        ids_pos <- (out$dim[id_end:id_start_new] - 1) * K + out$dim[ii]
        vtemp <- rep(0, length(y_evals))
        for (id_m in 1:L) {
          vtemp <- vtemp + dnsbeta(y_evals, as_array[i+1, id_m, ids_pos],
                                   bs_array[i+1, id_m, ids_pos], 0, T_kernel)*ps_array[i+1, id_m, ids_pos]*(1-epsilon_mcmc[i+1]) +
            dnsbeta(y_evals, a0_mcmc[i+1, id_m] , b0_mcmc[i+1, id_m], 0, T_kernel)*p0_mcmc[i+1, id_m]*(epsilon_mcmc[i+1])
        }
        ptemp[2:(ids_nonzero + 1)] <- ptemp[2:(ids_nonzero + 1)] + log(vtemp)
        ptemp <- ptemp[1:(ids_nonzero + 1)]
        res_temp <- exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp)))
        sampled <- sample(1:length(res_temp), size = 1, prob = res_temp)
        B[ii,] <- 0
        B[ii,sampled] <- 1#sample(1:min(col_max,ii), size = 1, prob = res_temp)
      }
    }

    # update Z again ?
    for (id_i in 1:K) {
      for (id_j in 1:K) {
        pos <- (id_i - 1) * K + id_j
        dists_B <- B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)]
        dists_y <- y[ids[[pos]]$Var2] - y[ids[[pos]]$Var1]
        interarr <- dists_B * dists_y
        y_DP <- interarr[interarr > 0]
        n_DP <- length(y_DP)
        for (j in 1:n_DP) {
          if (PRIOR_SPEC == "ab") {
            log_prob_vec1 <- log(p0_mcmc[i+1,]) + dnsbeta(x = y_DP[j],
                                                          shape1 = a0_mcmc[i+1,],
                                                          shape2 = b0_mcmc[i+1,],
                                                          min = 0, max = T_kernel, log = TRUE) +
              log(epsilon_mcmc[i])
            log_prob_vec2 <- log(p_mcmc[[(id_i-1)*K + id_j]][i+1,]) + dnsbeta(x = y_DP[j],
                                                                              shape1 = as_mcmc[[(id_i-1)*K + id_j]][i+1,],
                                                                              shape2 = bs_mcmc[[(id_i-1)*K + id_j]][i+1,],
                                                                              min = 0, max = T_kernel, log = TRUE) +
              log(1-epsilon_mcmc[i])
            log_prob_vec <- c(log_prob_vec1, log_prob_vec2)
          } else {
            log_prob_vec <- log(p_mcmc[i+1,]) + dnsbeta(x = y_DP[j],
                                                        shape1 = ss_mcmc[i,]*ms_mcmc[i,],
                                                        shape2 = ss_mcmc[i,]*(1-ms_mcmc[i+1,]),
                                                        min = 0, max = T_kernel, log = TRUE)
          }
          if (STORE_Z) {
            Z_mcmc[[(id_i-1)*K + id_j]][i+1,j] <- rcat_lpw(log_prob_vec)
          } else {
            Z_mcmc[[(id_i-1)*K + id_j]][j] <- rcat_lpw(log_prob_vec)
          }
        }
        if (STORE_Z) {
          Z_mcmc[[(id_i-1)*K + id_j]][i+1,n_DP:N] <- 0
        } else {
          Z_mcmc[[(id_i-1)*K + id_j]][(n_DP+1):N] <- 0
        }
      }
    }
    alphas_mcmc[[i+1]] <- alpha
    mus_mcmc[[i+1]] <- mu
    if (i %% 50 == 0) {
      print(i)
    }
  }
  fileout <- list(
    p_mcmc = p_mcmc, as_mcmc = as_mcmc,
    bs_mcmc = bs_mcmc, p0_mcmc = p0_mcmc,
    a0_mcmc = a0_mcmc, b0_mcmc = b0_mcmc,
    mus_mcmc = mus_mcmc, alphas_mcmc = alphas_mcmc,
    epsilon_mcmc = epsilon_mcmc,
    time = as.numeric(Sys.time() - start_time)
  )
  return(fileout)
}
#
