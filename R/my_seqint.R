# HELPER FUNCTION FOR GENERATE THE MATRIX
my_seqint <- function(x) {
  seq_from <- x[1]
  seq_to <- x[2]
  return(seq.int(from = seq_from, to = seq_to))
}
