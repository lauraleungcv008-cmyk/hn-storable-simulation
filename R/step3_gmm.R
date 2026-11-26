# Module 3(b): step 3 estimated by GMM instead of MLE
#
# Six moments of the size choice, each a ratio of household-level sums:
#   m1 Pr(buy | at least one brand on sale)     m2 Pr(buy | no sale)
#   m3 Pr(large | buy, sale)                    m4 Pr(large | buy, no sale)
#   m5 Pr(buy next week | bought large)         m6 Pr(buy next week | bought small)
# Model moments are computed exactly from the ergodic distribution of the
# solved model, so the objective is smooth. L = 6, K = 3: over-identified.
#
#   theta_hat = argmin g' W g,  g = m_data - m_model(theta),  W = Var(m_data)^{-1}

MOMENT_NAMES <- c("Pr(buy | sale)", "Pr(buy | no sale)",
                  "Pr(large | buy, sale)", "Pr(large | buy, no sale)",
                  "Pr(buy next | bought large)", "Pr(buy next | bought small)")

data_moments <- function(m, panel) {
  sale <- (m$states$sA[panel$K] + m$states$sB[panel$K]) > 0
  dim(sale) <- dim(panel$K)
  buy <- panel$D > 0
  large <- panel$D == 2
  small <- panel$D == 1
  T <- ncol(panel$K)
  b0 <- buy[, -T]; b1 <- buy[, -1]; l0 <- large[, -T]; s0 <- small[, -T]
  a <- cbind(rowSums(sale & buy), rowSums(!sale & buy),
             rowSums(sale & large), rowSums(!sale & large),
             rowSums(l0 & b1), rowSums(s0 & b1))
  b <- cbind(rowSums(sale), rowSums(!sale),
             rowSums(sale & buy), rowSums(!sale & buy),
             rowSums(l0), rowSums(s0))
  H <- nrow(a)
  A <- colSums(a); Bs <- colSums(b)
  mom <- A / Bs
  psi <- sweep(a - sweep(b, 2, mom, "*"), 2, Bs / H, "/")   # influence functions
  Sigma <- crossprod(psi) / H^2 * H / (H - 1)                # Var(m_hat), clustered
  list(m = mom, Sigma = Sigma)
}

model_moments <- function(m, theta3, omega, Pi) {
  P <- solve_size_problem(m, theta3, omega, Pi)
  mu <- stationary(m, P, Pi)
  sale <- (m$states$sA + m$states$sB) > 0
  buy <- 1 - P[, , 1]
  # Pr(buy next week | (k, i), size d): sum_k' Pi[k, k'] * buy[k', i'(i, d)]
  next_buy <- function(d) {
    nb <- Pi %*% buy[, m$i_next[, d + 1] + 1]       # nK x nI
    sum(mu * P[, , d + 1] * nb) / sum(mu * P[, , d + 1])
  }
  c(sum(mu[sale, ] * buy[sale, ]) / sum(mu[sale, ]),
    sum(mu[!sale, ] * buy[!sale, ]) / sum(mu[!sale, ]),
    sum(mu[sale, ] * P[sale, , 3]) / sum(mu[sale, ] * buy[sale, ]),
    sum(mu[!sale, ] * P[!sale, , 3]) / sum(mu[!sale, ] * buy[!sale, ]),
    next_buy(2), next_buy(1))
}

estimate_step3_gmm <- function(m, panel, omega, Pi, spec = "common", start = NULL) {
  dm <- data_moments(m, panel)
  W <- solve(dm$Sigma)
  K <- length(PAR_NAMES[[spec]])
  if (is.null(start)) start <- c(rep(0.2, K - 1), 0.1)
  unz <- function(z) c(z[-K], exp(z[K]))
  mom <- function(par) model_moments(m, to_theta3(par, spec), omega, Pi)
  obj <- function(par) {
    q <- tryCatch({ g <- dm$m - mom(par); drop(t(g) %*% W %*% g) }, error = function(e) NA)
    if (is.finite(q)) q else 1e12
  }
  fit <- optim(c(start[-K], log(start[K])), function(z) obj(unz(z)), method = "BFGS",
               control = list(reltol = 1e-14, maxit = 1000))
  par <- unz(fit$par)
  # Jacobian of the model moments, then the efficient-GMM variance (G' W G)^{-1}
  G <- sapply(1:K, function(j) {
    h <- 1e-5 * max(abs(par[j]), 1e-3)
    tp <- par; tm <- par; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (mom(tp) - mom(tm)) / (2 * h)
  })
  V <- tryCatch(solve(t(G) %*% W %*% G), error = function(e) matrix(NA, K, K))
  nm <- PAR_NAMES[[spec]]
  df <- 6 - K
  list(theta = setNames(par, nm), se = setNames(sqrt(pmax(diag(V), 0)), nm), G = G,
       J = fit$value, J_df = df, J_pvalue = 1 - pchisq(fit$value, df),
       converged = fit$convergence == 0)
}
