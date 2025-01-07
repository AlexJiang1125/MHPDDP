# FUNCTION: using the approximation trick in Ma et al. (2011) to calculate
# the expected log beta distribution over a and b
# returns E_{a,b}[Beta(a,b)], with inputs a,b
exp_lbeta <- function(ca,cb,da,db) {
  u_bar <- ca/da
  v_bar <- cb/db
  # E[ln(a)]
  qla <- digamma(ca) - log(da)
  # E[ln(b)]
  qlb <- digamma(cb) - log(db)
  # E[(ln(u) - ln(ubar))^2]
  qlusq <- (digamma(ca) - log(ca))^2 + psigamma(ca, deriv = 1)
  # E[(ln(v) - ln(vbar))^2]
  qlvsq <- (digamma(cb) - log(cb))^2 + psigamma(cb, deriv = 1)
  return(-lbeta(u_bar, v_bar) +
           u_bar*(digamma(u_bar + v_bar) - digamma(u_bar))*(qla - log(u_bar)) +
           v_bar*(digamma(u_bar + v_bar) - digamma(v_bar))*(qlb - log(v_bar)) +
           0.5*u_bar^2*(psigamma(u_bar + v_bar, deriv = 1) - psigamma(u_bar, deriv = 1))*qlusq +
           0.5*v_bar^2*(psigamma(u_bar + v_bar, deriv = 1) - psigamma(v_bar, deriv = 1))*qlvsq +
           u_bar*v_bar*psigamma(u_bar + v_bar, deriv = 1)*(qla - log(u_bar))*(qlb - log(v_bar)))
}
