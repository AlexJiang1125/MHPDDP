# load the parameters for Dirichlet processes
get_DPpars <- function(EXAMPLE_ID, samepriors, epsilon_true = 0.6, spec = "identity") {
  alpha_DP <- 1 # DP concentration parameter
  if (EXAMPLE_ID == 6) {
    T_kernel <- 5
  } else if (EXAMPLE_ID == 8) {
    T_kernel <- 5
  } else {
    T_kernel <- 1
  }

  # DP parameters
  if (EXAMPLE_ID == 1) {
    # common component
    p0s <- c(1)
    a0s <- c(1.2)
    b0s <- c(1.2)
    # idiosyncratic components
    ps <- list(
      list(c(0.4, 0.6),
           c(0.6, 0.4)),
      list(c(0.6, 0.4),
           c(0.4, 0.6))
    )

    as <- list(
      list(c(6, 3),
           c(3, 6)),
      list(c(4, 8),
           c(8, 4))
    )

    bs <- list(
      list(c(2, 8),
           c(8, 2)),
      list(c(3, 9),
           c(9, 3))
    )


  } else if (EXAMPLE_ID == 2) {
    # common component
    p0s <- c(1)
    a0s <- c(0.5)
    b0s <- c(2)

    # idiosyncratic components
    ps <- list(
      list(c(0.4, 0.6),
           c(0.6, 0.4)),
      list(c(0.6, 0.4),
           c(0.4, 0.6))
    )

    as <- list(
      list(c(3, 2),
           c(2, 8)),
      list(c(10, 2),
           c(4, 5))
    )

    bs <- list(
      list(c(12, 15),
           c(15, 8)),
      list(c(5, 10),
           c(10, 15))
    )
  } else if (EXAMPLE_ID == 3) {
    # common component
    p0s <- c(1)
    a0s <- c(0.5)
    b0s <- c(2)

    # idiosyncratic components
    ps <- list(
      list(c(0.4, 0.6),
           c(0.6, 0.4)),
      list(c(0.6, 0.4),
           c(0.4, 0.6))
    )

    as <- list(
      list(c(15, 12),
           c(12, 15)),
      list(c(15, 10),
           c(10, 15))
    )

    bs <- list(
      list(c(3, 5),
           c(5, 3)),
      list(c(3, 2),
           c(2, 6))
    )
  } else if (EXAMPLE_ID == 4) {
    p0s <- 1
    a0s <- 1
    b0s <- 4

    ps <- list(list(1,1), list(1,1))
    as <- list(list(2,4), list(6,1))
    bs <- list(list(6,1), list(2,1))
  } else if (EXAMPLE_ID == "flipped") {
    p0s <- 1
    a0s <- 1
    b0s <- 4

    ps <- list(list(1,1), list(1,1))
    as <- list(list(1,6), list(4,2))
    bs <- list(list(1,2), list(1,6))
  } else if (EXAMPLE_ID == 5) {
    # triangle distribution
    # mixture weights between the two triangles
    eps_tri0 <- eps_tri <- 0.9
    # tri_abc refers to min, mode and max
    ## common parts
    tri_a0 <- c(0, 0.4)
    tri_b0 <- c(0.02, 0.52)
    tri_c0 <- c(0.5, 1)
    ## idio parts
    # the (2,2) par is Beta distribution!
    tri_a <- list(c(0, 0.3),
                  c(0.5, 0),
                  c(0.4,0),
                  c(1,1))
    tri_b <- list(c(0.15, 0.52),
                  c(0.98, 0.48),
                  c(0.85,0.48),
                  c(1,1))
    tri_c <- list(c(0.6, 1),
                  c(1, 0.6),
                  c(1,0.7),
                  c(1,1))

    p0s <- 1
    a0s <- 1
    b0s <- 4

    ps <- list(list(1,1), list(1,1))
    as <- list(list(2,4), list(6,1))
    bs <- list(list(6,1), list(2,1))
  } else if (EXAMPLE_ID == 6) {
    exp_0 <- 1
    ps <- list(list(1,1), list(1,1))
    exp_s <- list(list(0.5, 2), list(2, 0.5))

    p0s <- 1
    a0s <- 1
    b0s <- 4

    ps <- list(list(1,1), list(1,1))
    as <- list(list(2,4), list(6,1))
    bs <- list(list(6,1), list(2,1))
  } else if (EXAMPLE_ID == 7) {
    load("Data/pars/BetaDPM_data_common.RData")
    list2env(res, envir = globalenv())
    p0s <- dir_weights
    a0s <- a_dir
    b0s <- b_dir
    load("Data/pars/BetaDPM_data_11.RData")
    res11 <- res
    load("Data/pars/BetaDPM_data_12.RData")
    res12 <- res
    load("Data/pars/BetaDPM_data_21.RData")
    res21 <- res
    load("Data/pars/BetaDPM_data_22.RData")
    res22 <- res
    ps <- list(list(res11$dir_weights, res12$dir_weights),
               list(res21$dir_weights, res22$dir_weights))
    as <- list(list(res11$a_dir, res12$a_dir),
               list(res21$a_dir, res22$a_dir))
    bs <- list(list(res11$b_dir, res12$b_dir),
               list(res21$b_dir, res22$b_dir))
  } else if (EXAMPLE_ID == 8) {
    exp_0 <- 1
    ps <- list(list(1,1), list(1,1))
    exp_s <- list(list(2, 0.8), list(0.8, 2))

    p0s <- 1
    a0s <- 1
    b0s <- 4

    ps <- list(list(1,1), list(1,1))
    as <- list(list(2,4), list(6,1))
    bs <- list(list(6,1), list(2,1))
  } else if (EXAMPLE_ID == 9) {
    p0s <- 1
    a0s <- 1
    b0s <- 4

    ps <- list(list(1,1), list(1,1))
    as <- list(list(2,4), list(1.5,1))
    bs <- list(list(6,1), list(5,1))
  }

  # set ca/cb based on choice of a,b values

  # prior hyperparameters
  if (samepriors) {
    a_mean <- mean(c(unlist(as), a0s))
    b_mean <- mean(c(unlist(bs), b0s))
    d_a <- d_b <- d_a0 <- d_b0 <- 1
    if (spec == "identity") {
      c_a0 <- c_a <- a_mean
      c_b0 <- c_b <- b_mean
    } else {
      c_a0 <- c_a <- max(0.5*(-1 + sqrt(1 + 4*a_mean)),1)
      c_b0 <- c_b <- 0.5*(-1 + sqrt(1 + 4*b_mean))
    }
  } else {
    a_mean <- mean(unlist(as))
    b_mean <- mean(unlist(bs))
    c_a <- a_mean
    c_b <- b_mean
    d_a <- 1
    d_b <- 1
    # check
    # hist(rgamma(1000, c_b, d_b))
    # abline(v = b_mean)
    c_a0 <- c_b0 <- 2
    d_a0 <- d_b0 <- 2
    # hist(rgamma(1000, c_a0, d_a0))
    # abline(v = a0s)
  }

  if (EXAMPLE_ID == 5) {
    return(
      res = list(
        ps = ps, as = as, bs = bs,
        p0s = p0s, a0s = a0s, b0s = b0s,
        c_a = c_a, c_b = c_b,
        c_a0 = c_a0, c_b0 = c_b0,
        d_a = d_a, d_b = d_b,
        d_a0 = d_a0, d_b0 = d_b0,
        alpha_DP = alpha_DP,
        epsilon_true = epsilon_true,
        eps_tri = eps_tri, eps_tri0 = eps_tri0,
        tri_a0 = tri_a0, tri_b0 = tri_b0, tri_c0 = tri_c0,
        tri_a = tri_a, tri_b = tri_b, tri_c = tri_c,
        T_kernel = T_kernel
      )
    )
  } else if (EXAMPLE_ID %in% c(6,8)) {
    return(
      res = list(
        ps = ps, as = as, bs = bs,
        p0s = p0s, a0s = a0s, b0s = b0s,
        c_a = c_a, c_b = c_b,
        c_a0 = c_a0, c_b0 = c_b0,
        d_a = d_a, d_b = d_b,
        d_a0 = d_a0, d_b0 = d_b0,
        alpha_DP = alpha_DP,
        epsilon_true = epsilon_true,
        T_kernel = T_kernel,
        exp_0 = exp_0, ps = ps, exp_s = exp_s
      )
    )
  } else {
    return(
      res = list(
        ps = ps, as = as, bs = bs,
        p0s = p0s, a0s = a0s, b0s = b0s,
        c_a = c_a, c_b = c_b,
        c_a0 = c_a0, c_b0 = c_b0,
        d_a = d_a, d_b = d_b,
        d_a0 = d_a0, d_b0 = d_b0,
        alpha_DP = alpha_DP,
        epsilon_true = epsilon_true,
        T_kernel = T_kernel
      )
    )
  }
}

