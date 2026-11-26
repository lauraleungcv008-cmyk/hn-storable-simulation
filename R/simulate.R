# Module 1 (continued): simulate a household purchase panel
#
# The estimators see only: price state K[h, t] (hence all prices), size
# choice D[h, t] (0 none, 1 small, 2 large) and brand B[h, t] (0 none, 1 A, 2 B).
# Inventory INV[h, t] is returned for diagnostics only and is never used
# in estimation, as in real scanner data.

simulate_panel <- function(m, P, alpha, xi_B, H, T, seed) {
  set.seed(seed)
  mu <- stationary(m, P)
  s0 <- sample.int(m$nK * m$nI, H, replace = TRUE, prob = as.vector(mu))
  k <- (s0 - 1) %% m$nK + 1
  i <- (s0 - 1) %/% m$nK                                   # 0..I_max

  # Pr(brand B | size s, price state k): static logit within size
  vA <- -alpha * m$prices[, 1, ]
  vB <- xi_B - alpha * m$prices[, 2, ]
  pB <- exp(vB) / (exp(vA) + exp(vB))                      # nK x 2

  K <- D <- B <- INV <- matrix(0L, H, T)
  for (t in 1:T) {
    K[, t] <- k
    INV[, t] <- i
    p0 <- P[cbind(k, i + 1, 1)]
    p1 <- P[cbind(k, i + 1, 2)]
    draw <- runif(H)
    d <- as.integer(draw > p0) + as.integer(draw > p0 + p1)  # 0, 1, 2
    D[, t] <- d
    brand <- integer(H)
    buy <- d > 0
    brand[buy] <- 1L + as.integer(runif(sum(buy)) < pB[cbind(k[buy], d[buy])])
    B[, t] <- brand
    i <- m$i_next[cbind(i + 1, d + 1)]
    sA <- as.integer(runif(H) < m$Pi1[m$states$sA[k] + 1, 2])
    sB <- as.integer(runif(H) < m$Pi1[m$states$sB[k] + 1, 2])
    k <- sA + 2L * sB + 1L
  }
  list(K = K, D = D, B = B, INV = INV)
}
