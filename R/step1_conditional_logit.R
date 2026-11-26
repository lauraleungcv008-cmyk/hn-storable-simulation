# Step 1 of Hendel and Nevo (2006): brand choice conditional on size
#
# On each purchase occasion, restrict the choice set to the brands of the
# size actually bought and estimate a static conditional logit by maximum
# likelihood. No dynamic programming is needed, because brand choice does
# not change future inventory.
#
#   Pr(brand b | size s, week t) = exp(xi_b - alpha * p_bst) / sum_b' exp(xi_b' - alpha * p_b'st)
#
# Parameters: theta1 = (alpha, xi_B), with xi_A = 0.

step1_data <- function(m, panel) {
  buy <- panel$D > 0
  k <- panel$K[buy]
  s <- panel$D[buy]
  list(pA = m$prices[cbind(k, 1, s)],
       pB = m$prices[cbind(k, 2, s)],
       chose_B = as.numeric(panel$B[buy] == 2))
}

# Log-likelihood and its gradient
step1_loglik <- function(theta1, dat) {
  alpha <- theta1[1]; xi_B <- theta1[2]
  vA <- -alpha * dat$pA
  vB <- xi_B - alpha * dat$pB
  mx <- pmax(vA, vB)
  log_denom <- mx + log(exp(vA - mx) + exp(vB - mx))
  ll <- sum(dat$chose_B * vB + (1 - dat$chose_B) * vA - log_denom)
  probB <- exp(vB - log_denom)
  resid <- dat$chose_B - probB
  grad <- c(sum(-(dat$pB - dat$pA) * resid),    # d ll / d alpha
            sum(resid))                          # d ll / d xi_B
  attr(ll, "gradient") <- grad
  ll
}

estimate_step1 <- function(m, panel, start = c(1, 0)) {
  dat <- step1_data(m, panel)
  fn <- function(th) -step1_loglik(th, dat)
  gr <- function(th) -attr(step1_loglik(th, dat), "gradient")
  fit <- optim(start, fn, gr, method = "BFGS", hessian = TRUE,
               control = list(reltol = 1e-12, maxit = 500))
  V <- solve(fit$hessian)
  list(theta = setNames(fit$par, c("alpha", "xi_B")),
       se = setNames(sqrt(diag(V)), c("alpha", "xi_B")),
       V = V, n_purchases = length(dat$chose_B), converged = fit$convergence == 0)
}
