listRows<-function(m){
  #converts a sparse Matrix into a list of its columns
  m = t(m)
  res<-split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
  names(res)<-rownames(m)
  return(res)
}
