exp_norm_rmzero <- function(ptemp){
  id_nonzero <- which(ptemp != 0)
  ptemp_nonzero <- ptemp[id_nonzero]
  res_nonzero <- exp(ptemp_nonzero - max(ptemp_nonzero))/sum(exp(ptemp_nonzero - max(ptemp_nonzero)))
  res <- rep(0 , length(ptemp))
  res[id_nonzero] <- res_nonzero
  return(res)
}
