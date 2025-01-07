# FUNCTION: sample from categorical distribution based on the weights
rcat_lpw <- function(p) {
  lp <- p - max(p)
  p_normalized <- exp(lp)/(sum(exp(lp)))
  return(sample(1:length(p_normalized), size = 1, prob = p_normalized))
}
