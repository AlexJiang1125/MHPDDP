# FUNCTION: helper function for exponentiated normalization
exp_norm <- function(ptemp){
  return(exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp))))
}
