#' @import MCMCpack
#' @import extraDistr
#' @export
MHPDDP_svi <- function(data = df,
                       K = 2, L = 15,
                       split_ratio = 0.1,
                       e_alpha = 2, f_alpha = 4,
                       T_all = 15000, a_mu = 1, b_mu = 100,
                       T_kernel = 1, c_a = 2, c_b = 3.5,
                       d_a = 1, d_b = 1, alpha_DP = 1,
                       c_a0 = 2, c_b0 = 2, d_a0 = 2, d_b0 = 2,
                       seed = 123, borrow_mode = "RANDOM",
                       sgd_par1 = 1, sgd_par2 = 0.51
                       ) {



  # read in the data
  y <- data$timestamp
  N <- length(y)

  start_time <- Sys.time()

  col_max = 200
  # get transition pairs
  ids_all <- list()
  pos <- 1
  for (j in 1:K) {
    for (k in 1:K) {
      id_j <- which(data$dim == j)
      id_k <- which(data$dim == k)
      id_jk <- expand.grid(id_j, id_k)
      id_jk <- id_jk[id_jk$Var2 - id_jk$Var1 > 0, ]
      id_jk$time_start <- data$timestamp[id_jk$Var1]
      id_jk$time_end <- data$timestamp[id_jk$Var2]
      #id_jk <- id_jk[out_all$parentdim[id_jk$Var2] != 0 , ]
      id_jk <- id_jk[id_jk$time_end - id_jk$time_start < T_kernel, ]
      ids_all[[pos]] <- id_jk
      pos <- pos + 1
    }
  }

  # variational parameters for common components
  ca0_vi <- rep(c_a0, L) + runif(L, -1, 1) # rgamma(L, 40, 20)  #matrix(0, nrow = B, ncol = L)
  cb0_vi <- rep(c_b0, L) + runif(L, -1, 1) # rgamma(L, 40, 20)  #matrix(0, nrow = B, ncol = L)
  da0_vi <- rep(d_a0, L) + runif(L, -0.2, 0.2) # rgamma(L, 10, 50)  #matrix(0, nrow = B, ncol = L)
  db0_vi <- rep(d_b0, L) + runif(L, -0.2, 0.2) # rgamma(L, 10, 50)  #matrix(0, nrow = B, ncol = L)

  # variational parameters for idiosyncratic components
  # like MCMC alg we use a list of vectors to store
  # for dimension (id_i) to transit to (id_j), it is [[(id_i-1)*K+id_j]]
  cas_vi <- list()
  cbs_vi <- list()
  das_vi <- list()
  dbs_vi <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      cas_vi[[(id_i - 1)*K + id_j]] <- rep(c_a, L) + runif(L, -1, 1) #rgamma(L, 40, 20)
      cbs_vi[[(id_i - 1)*K + id_j]] <- rep(c_b, L) + runif(L, -1, 1) #rgamma(L, 40, 20)
      das_vi[[(id_i - 1)*K + id_j]] <- rep(d_a, L) + runif(L, -0.2, 0.2) #rgamma(L, 10, 50)
      dbs_vi[[(id_i - 1)*K + id_j]] <- rep(d_b, L) + runif(L, -0.2, 0.2) #rgamma(L, 10, 50)

    }
  }

  ca0_prior <- rep(c_a, L)
  cb0_prior <- rep(c_b, L)
  da0_prior <- rep(d_a, L)
  db0_prior <- rep(d_b, L)
  cas_prior <- cbs_prior <- das_prior <- dbs_prior <- list()
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      cas_prior[[(id_i - 1)*K + id_j]] <- rep(c_a, L) #rgamma(L, 40, 20)
      cbs_prior[[(id_i - 1)*K + id_j]] <- rep(c_b, L) #rgamma(L, 40, 20)
      das_prior[[(id_i - 1)*K + id_j]] <- rep(d_a, L) #rgamma(L, 10, 50)
      dbs_prior[[(id_i - 1)*K + id_j]] <- rep(d_b, L) #rgamma(L, 10, 50)
    }
  }
  # variational parameters for p_l
  alpha0_vi <- rgamma(L, alpha_DP, L) # rep(alpha_DP/L, L)  #matrix(0, nrow = B, ncol = L)
  alphas_vi <- alphas_vi_init <- list()
  Zs_vi <- list()

  alpha0_vi_init <- rep(alpha_DP/L, L)

  # variational parameters for allocation indicators
  for (id_i in 1:K) {
    for (id_j in 1:K) {
      alphas_vi[[(id_i - 1)*K + id_j]] <- rgamma(L, alpha_DP, L) # rep(alpha_DP/L, L)
      alphas_vi_init[[(id_i - 1)*K + id_j]] <- rep(alpha_DP/L, L)
    }
  }

  eps0_vi <- c(1, 1)
  tp <- runif(1, 0, 1)
  eps_vi <- c(tp, 1-tp)
  ### HP parameters ###
  # dimensional transition from dim id_i to id_j is stored
  # in the (id_i - 1)*K + id_j
  eta_alphas1 <- eta_alphas2 <- list()
  # eta_alphas1 has mean 2
  eta_alphas1_prior <- eta_alphas1 <- rep(e_alpha, K^2)  #rgamma(K^2, 4, 2)
  # eta_alphas2 has mean 4
  eta_alphas2_prior <- eta_alphas2 <- rep(f_alpha, K^2)  #rgamma(K^2, 8, 2)
  # eta_mus1 has mean 1
  eta_mus1_prior <- eta_mus1 <- rep(a_mu,K) #rgamma(K, 2, 2)
  # eta_mus1 has mean 4
  eta_mus2_prior <- eta_mus2 <- rep(b_mu,K) #rgamma(K, 15, 1)

  eta_alphas1_init <- eta_alphas1 + runif(K^2, min = -1, max = 1)
  eta_alphas2_init <- eta_alphas2 + runif(K^2, min = -2, max = 2)
  eta_mus1_init <- eta_mus1 + runif(K, min = -0.2, max = 0.2)
  eta_mus2_init <- eta_mus2 + runif(K, min = -5, max = 5)

  # eta_alphas1_init <- as.vector(alpha_true)*1000
  # eta_alphas2_init <- c(100, 100, 100, 100)

  # rel.err
  rel_err <- 1
  us_bar <- vs_bar <- list()
  qlpis <- qlas <- qlbs <- qlusqs <- qlvsqs <- list()

  ######## VI ALGORITHM #######
  start_time = Sys.time()

  i <- 1
  rmises <- NA
  elbos <- NA
  elbo_pos <- 1
  aocs <- NA
  pos_i <- 1
  RUN <- TRUE

  #############################

  #while (RUN && as.numeric(difftime(Sys.time(),start_time,units="secs")) < 1800) {
  #
  #}
  while (as.numeric(difftime(Sys.time(),start_time,units="secs")) < 10) {
    i <- i + 1
    rho <- 1*(i + sgd_par1) ^ (-sgd_par2)
    # GENERATE SUBSAMPLES
    res <- generate_subsample(data, ids_all, split_ratio, T_all, K)
    out <- res$out
    ids <- res$ids
    N <- dim(out)[1]
    Ns <- sapply(ids, FUN = "dim")[1,]
    Ns_alpha <- as.vector(table(out$dim))
    B <- matrix(0,N,N)
    diag(B) <- 1

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
            eps_part <- 1*(digamma(eps_vi[1]) - digamma(eps_vi[1] + eps_vi[2]))
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
            eps_part <- 1*(digamma(eps_vi[2]) - digamma(eps_vi[1] + eps_vi[2]))
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
  }
}
