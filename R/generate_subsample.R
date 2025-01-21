# FUNCTI_CLUS: generate subsample from the whole dataset
# the input is a list of vectors, each having N
# in order to scale with the HP version which we do next
# here we sample a 'random window' and then select elements
# that falls within the random window
# we also split the info matrix about indexs
generate_subsample <- function(out_all, ids_all, split_ratio, T_all, K) {
  # generate random start and end dates
  # T_all <- length(y[[1]])
  time_start <- runif(1, 0, T_all - T_all * split_ratio)
  time_end <- time_start + T_all * split_ratio
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
  res <- list(
    out = out,
    ids = ids
  )
  return(res)
}
