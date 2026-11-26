# Step 2 of Hendel and Nevo (2006): inclusive values and their transitions
#
# omega_s(k) summarises all brands of size s into one quality-adjusted price
# index. With few discrete price states, omega takes one value per state k,
# so its first-order Markov transition is estimated by transition frequencies.
# (HN have many brands and continuous prices; they regress omega_t on omega_{t-1}.)

estimate_step2 <- function(m, panel, theta1) {
  omega <- inclusive_values(m, theta1["alpha"], theta1["xi_B"])
  K <- panel$K
  from <- as.vector(K[, -ncol(K)])
  to <- as.vector(K[, -1])
  counts <- table(factor(from, 1:m$nK), factor(to, 1:m$nK))
  Pi_hat <- unclass(counts) / rowSums(counts)
  list(omega = omega, Pi = matrix(Pi_hat, m$nK, m$nK))
}
