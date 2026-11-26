# Module 1: a stripped-down Hendel-Nevo (2006) storable-good model
#
# Products: 2 brands (A, B) x 2 sizes (small = 2 units, large = 4 units).
# Prices:   each brand is on regular price or on sale; the two brands follow
#           independent two-state Markov chains. Joint price state k = 1..4.
# Household: consumes 1 unit per week if it has one; holds 0..I_max units.
#           Each week it picks one of 5 options: no purchase, or (brand, size).
#
# Utility of buying brand b in size s in week t:
#     xi_b + eta_s - alpha * p_bst - c * i' + eps_bs
# Utility of no purchase:
#     - c * i' - kappa * 1[stock-out] + eps_0
# i' = min(max(i + q_s - 1, 0), I_max); eps are i.i.d. Type I extreme value.
#
# HN's key restrictions hold by construction: brand matters only at purchase
# and storage counts units regardless of brand. Hence the 5-option problem
# collapses to a 3-option size problem in which each size carries its
# inclusive value omega_s(k) = log sum_b exp(xi_b - alpha * p_bs(k)).

EULER <- 0.5772156649015329

make_model <- function(beta = 0.95, I_max = 12, q = c(2, 4),
                       unit_price = c(1.0, 0.7), size_discount = c(1.0, 0.9),
                       Pi1 = matrix(c(0.75, 0.25, 0.70, 0.30), 2, byrow = TRUE),
                       kappa = 6) {
  states <- expand.grid(sA = 0:1, sB = 0:1)          # row k: 0 = regular, 1 = sale
  nK <- nrow(states)
  Pi <- matrix(0, nK, nK)                            # Pi[k, k'] = Pr(k' | k)
  for (k in 1:nK) for (kk in 1:nK)
    Pi[k, kk] <- Pi1[states$sA[k] + 1, states$sA[kk] + 1] *
                 Pi1[states$sB[k] + 1, states$sB[kk] + 1]
  prices <- array(0, c(nK, 2, 2))                    # prices[k, brand, size]
  for (k in 1:nK) for (b in 1:2) for (s in 1:2) {
    on_sale <- c(states$sA[k], states$sB[k])[b]
    prices[k, b, s] <- q[s] * unit_price[on_sale + 1] * size_discount[s]
  }
  inv <- 0:I_max
  qd <- c(0, q)                                      # units bought by size choice d = 0, 1, 2
  i_next   <- sapply(qd, function(x) pmin(pmax(inv + x - 1, 0), I_max))   # nI x 3
  stockout <- sapply(qd, function(x) as.numeric(inv + x == 0))            # nI x 3
  list(beta = beta, I_max = I_max, nI = I_max + 1, nK = nK, q = q, qd = qd,
       Pi = Pi, Pi1 = Pi1, prices = prices, states = states,
       i_next = i_next, stockout = stockout, kappa = kappa)
}

# omega[k, s]: inclusive value of size s in price state k
inclusive_values <- function(m, alpha, xi_B) {
  xi <- c(0, xi_B)                                   # brand A normalised to 0
  sapply(1:2, function(s)
    log(exp(xi[1] - alpha * m$prices[, 1, s]) + exp(xi[2] - alpha * m$prices[, 2, s])))
}

# Flow utility of the size problem, array [k, i, d], d = none / small / large
flow_utility <- function(m, omega, eta, c_store) {
  u <- array(0, c(m$nK, m$nI, 3))
  for (d in 1:3) {
    base <- -c_store * m$i_next[, d] - m$kappa * m$stockout[, d]          # length nI
    u[, , d] <- matrix(base, m$nK, m$nI, byrow = TRUE)
    if (d > 1) u[, , d] <- u[, , d] + eta[d - 1] + omega[, d - 1]
  }
  u
}

choice_values <- function(m, u, V, Pi = m$Pi) {
  EV <- Pi %*% V                                     # EV[k, i'] = E[V(k', i') | k]
  v <- u
  for (d in 1:3) v[, , d] <- u[, , d] + m$beta * EV[, m$i_next[, d] + 1]
  v
}

logsumexp3 <- function(v) {
  mx <- pmax(v[, , 1], v[, , 2], v[, , 3])
  mx + log(exp(v[, , 1] - mx) + exp(v[, , 2] - mx) + exp(v[, , 3] - mx))
}

softmax3 <- function(v) {
  lse <- logsumexp3(v)
  P <- v
  for (d in 1:3) P[, , d] <- exp(v[, , d] - lse)
  P
}

# Solution method 1: value function iteration (contraction with modulus beta)
solve_vfi <- function(m, u, Pi = m$Pi, tol = 1e-12, max_iter = 20000) {
  V <- matrix(0, m$nK, m$nI)
  for (it in 1:max_iter) {
    V_new <- logsumexp3(choice_values(m, u, V, Pi)) + EULER
    if (max(abs(V_new - V)) < tol) { V <- V_new; break }
    V <- V_new
  }
  list(V = V, P = softmax3(choice_values(m, u, V, Pi)), iter = it)
}

# Transition matrix over states s = (k, i), column-major: s = k + nK * i
transition_matrix <- function(m, P, Pi = m$Pi) {
  nS <- m$nK * m$nI
  F <- matrix(0, nS, nS)
  for (d in 1:3) for (i in 0:m$I_max) {
    rows <- (1:m$nK) + m$nK * i
    cols <- (1:m$nK) + m$nK * m$i_next[i + 1, d]
    F[rows, cols] <- F[rows, cols] + P[, i + 1, d] * Pi
  }
  F
}

# Solution method 2: policy iteration (Newton-Kantorovich steps)
solve_pi <- function(m, u, Pi = m$Pi, tol = 1e-13, max_iter = 200) {
  P <- array(1 / 3, dim(u))
  I_S <- diag(m$nK * m$nI)
  for (it in 1:max_iter) {
    e <- as.vector(rowSums(P * (u + EULER - log(pmax(P, 1e-300))), dims = 2))
    V <- matrix(solve(I_S - m$beta * transition_matrix(m, P, Pi), e), m$nK, m$nI)
    P_new <- softmax3(choice_values(m, u, V, Pi))
    if (max(abs(P_new - P)) < tol) { P <- P_new; break }
    P <- P_new
  }
  list(V = V, P = P, iter = it)
}

# Ergodic distribution over (k, i), returned as an nK x nI matrix
stationary <- function(m, P, Pi = m$Pi) {
  nS <- m$nK * m$nI
  A <- t(transition_matrix(m, P, Pi)) - diag(nS)
  A[nS, ] <- 1
  b <- c(rep(0, nS - 1), 1)
  matrix(solve(A, b), m$nK, m$nI)
}
