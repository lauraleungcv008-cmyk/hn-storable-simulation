# Step 3 of Hendel and Nevo (2006): dynamic size choice by nested fixed point MLE
#
# Given omega and its transition matrix from steps 1-2, the household solves a
# 3-option dynamic problem (none / small / large) with state (k, i).
# Parameters: theta3 = (eta_small, eta_large, c). c is optimised in logs.
#
# Unobserved inventory: weekly consumption is 1 unit, so the whole inventory
# path of a household is pinned down by its initial stock i0 and its observed
# purchases. The likelihood integrates over i0 using a prior distribution:
#   "stationary": model-implied ergodic distribution of i given the first-week
#                 price state (correct specification)
#   "empty":      every household starts with zero stock (misspecified)

# Count matrix N[(h, i0), cell]: number of weeks in which household h, if it
# started with stock i0, would be in state (k, i) and choose d.
# cell = k + nK * i + nK * nI * d. It does not depend on theta, so it is built
# once per data set.
build_counts <- function(m, panel) {
  H <- nrow(panel$K); T <- ncol(panel$K); nI <- m$nI
  nCell <- m$nK * nI * 3
  N <- matrix(0, H * nI, nCell)
  inv <- rep(0:m$I_max, each = H)                  # rows ordered h fastest, then i0
  rows <- seq_len(H * nI)
  for (t in 1:T) {
    k <- rep(panel$K[, t], nI)
    d <- rep(panel$D[, t], nI)
    cell <- k + m$nK * inv + m$nK * nI * d
    N[cbind(rows, cell)] <- N[cbind(rows, cell)] + 1
    inv <- m$i_next[cbind(inv + 1, d + 1)]
  }
  list(N = N, H = H, k1 = panel$K[, 1])
}

solve_size_problem <- function(m, theta3, omega, Pi) {
  u <- flow_utility(m, omega, eta = theta3[1:2], c_store = theta3[3])
  solve_pi(m, u, Pi)$P
}

step3_loglik <- function(theta3, m, counts, omega, Pi, prior = "stationary") {
  P <- solve_size_problem(m, theta3, omega, Pi)
  LL <- matrix(counts$N %*% log(as.vector(P)), counts$H, m$nI)   # H x nI (by i0)
  if (prior == "stationary") {
    mu <- stationary(m, P, Pi)                       # nK x nI
    w <- pmax(mu[counts$k1, , drop = FALSE], 0)
    w <- w / rowSums(w)
    mx <- apply(LL, 1, max)
    return(sum(mx + log(rowSums(w * exp(LL - mx)))))
  }
  sum(LL[, 1])                                       # "empty": i0 = 0 for everyone
}

# Parameter maps. "common": one purchase-occasion utility for both sizes,
# theta3 = (eta, eta, c). "free": separate eta_small and eta_large.
PAR_NAMES <- list(common = c("eta", "c"), free = c("eta_small", "eta_large", "c"))
to_theta3 <- function(par, spec)
  if (spec == "common") c(par[1], par[1], par[2]) else par

estimate_step3_mle <- function(m, counts, omega, Pi, spec = "common",
                               prior = "stationary", starts = NULL) {
  K <- length(PAR_NAMES[[spec]])
  # The likelihood is flat as c -> 0, and a single quasi-Newton run can drift
  # there and stop. Use bounded optimisation from two starting values and keep
  # the run with the highest likelihood.
  if (is.null(starts)) starts <- list(c(rep(0.3, K - 1), 0.03), c(rep(0.0, K - 1), 0.10))
  unz <- function(z) c(z[-K], exp(z[K]))               # c is optimised in logs
  nll <- function(par) {
    # Extreme trial values can make the inner solution fail; penalise them
    ll <- tryCatch(step3_loglik(to_theta3(par, spec), m, counts, omega, Pi, prior),
                   error = function(e) NA)
    if (is.finite(ll)) -ll else 1e12
  }
  lower <- c(rep(-5, K - 1), log(1e-3)); upper <- c(rep(5, K - 1), log(1))
  runs <- lapply(starts, function(s)
    optim(c(s[-K], log(s[K])), function(z) nll(unz(z)), method = "L-BFGS-B",
          lower = lower, upper = upper, control = list(factr = 1e5, maxit = 500)))
  fit <- runs[[which.min(sapply(runs, `[[`, "value"))]]
  par <- unz(fit$par)
  # Standard errors from the numerical Hessian. They treat omega and Pi as
  # known (no Murphy-Topel correction for steps 1-2).
  V <- tryCatch(solve(optimHess(par, nll)), error = function(e) matrix(NA, K, K))
  nm <- PAR_NAMES[[spec]]
  list(theta = setNames(par, nm), se = setNames(sqrt(pmax(diag(V), 0)), nm), V = V,
       loglik = -fit$value, converged = fit$convergence == 0,
       n_optima = length(unique(round(sapply(runs, `[[`, "value"), 2))))
}
