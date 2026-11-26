# Monte Carlo study of the three-step estimator
#   Rscript run_montecarlo.R [R_BASE]       [R_COMP]   (replications: E1-E2 default 200, E4 default 100)
#
# E1  Baseline: bias, RMSE and 95% CI coverage under correct specification;
#     step 3 by nested fixed point MLE and by GMM.
# E2  Restriction on purchase utility: step 3 with size-specific utilities
#     left free, on the same data.
# E4  Unobserved inventory and panel composition: N x T fixed at 208,000
#     household-weeks; correct versus wrong prior for initial inventory.
for (f in c("model", "simulate", "step1_conditional_logit", "step2_inclusive_values",
            "step3_mle", "step3_gmm", "truth")) source(file.path("R", paste0(f, ".R")))
library(parallel)
args <- as.integer(commandArgs(TRUE))
R_BASE <- if (length(args) >= 1) args[1] else 200     # E1, E2
R_COMP <- if (length(args) >= 2) args[2] else 100     # E4
CORES <- 2
m <- make_model()
omega0 <- inclusive_values(m, TRUTH["alpha"], TRUTH["xi_B"])
P0 <- solve_size_problem(m, TRUE_THETA3, omega0, m$Pi)

one_rep_baseline <- function(seed, H, T) {
  panel <- simulate_panel(m, P0, TRUTH["alpha"], TRUTH["xi_B"], H, T, seed)
  s1 <- estimate_step1(m, panel)
  s2 <- estimate_step2(m, panel, s1$theta)
  counts <- build_counts(m, panel)
  mle    <- estimate_step3_mle(m, counts, s2$omega, s2$Pi, "common")
  oracle <- estimate_step3_mle(m, counts, omega0, m$Pi, "common")   # true steps 1-2
  gmm    <- estimate_step3_gmm(m, panel, s2$omega, s2$Pi, "common")
  free   <- estimate_step3_mle(m, counts, s2$omega, s2$Pi, "free")
  c(s1 = s1$theta, s1se = s1$se, mle = mle$theta, mlese = mle$se,
    orc = oracle$theta, orcse = oracle$se, gmm = gmm$theta, gmmse = gmm$se,
    Jp = gmm$J_pvalue, free = free$theta,
    contrast = unname(free$theta["eta_large"] - 2 * free$theta["eta_small"]))
}

one_rep_composition <- function(seed, H, T) {
  panel <- simulate_panel(m, P0, TRUTH["alpha"], TRUTH["xi_B"], H, T, seed)
  s1 <- estimate_step1(m, panel)
  s2 <- estimate_step2(m, panel, s1$theta)
  counts <- build_counts(m, panel)
  a <- estimate_step3_mle(m, counts, s2$omega, s2$Pi, "common", prior = "stationary")
  b <- estimate_step3_mle(m, counts, s2$omega, s2$Pi, "common", prior = "empty")
  c(stat = a$theta, statse = a$se, empty = b$theta, emptyse = b$se)
}

run <- function(fun, H, T, tag, R = R_BASE) {
  seeds <- 1000 * H + T + seq_len(R)
  t0 <- Sys.time()
  out <- mclapply(seeds, fun, H = H, T = T, mc.cores = CORES)
  ok <- sapply(out, is.numeric)
  res <- do.call(rbind, out[ok])
  saveRDS(res, sprintf("results/%s_H%d_T%d.rds", tag, H, T))
  cat(sprintf("%s H=%d T=%d: %d/%d replications, %.0fs\n", tag, H, T, sum(ok), R,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res
}

for (H in c(500, 2000)) run(one_rep_baseline, H, 104, "baseline")
for (d in list(c(8000, 26), c(4000, 52), c(2000, 104), c(1000, 208)))
  run(one_rep_composition, d[1], d[2], "composition", R_COMP)
