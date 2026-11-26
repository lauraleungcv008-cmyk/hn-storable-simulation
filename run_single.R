# One simulated data set: numerical checks and the three-step estimator
#   Rscript run_single.R
for (f in c("model", "simulate", "step1_conditional_logit", "step2_inclusive_values",
            "step3_mle", "step3_gmm", "truth")) source(file.path("R", paste0(f, ".R")))

m <- make_model()
omega0 <- inclusive_values(m, TRUTH["alpha"], TRUTH["xi_B"])
u0 <- flow_utility(m, omega0, TRUE_THETA3[1:2], TRUE_THETA3[3])

# 1. Two solution methods agree
vfi <- solve_vfi(m, u0); pit <- solve_pi(m, u0); P0 <- pit$P
cat(sprintf("VFI: %d iterations | policy iteration: %d iterations | max |V_vfi - V_pi| = %.1e\n",
            vfi$iter, pit$iter, max(abs(vfi$V - pit$V))))

# 2. Simulated moments match exact model moments
big <- simulate_panel(m, P0, TRUTH["alpha"], TRUTH["xi_B"], H = 20000, T = 104, seed = 1)
mm <- model_moments(m, TRUE_THETA3, omega0, m$Pi); dm <- data_moments(m, big)$m
cat("\nMoment                          model   simulated (H = 20,000)\n")
for (j in 1:6) cat(sprintf("%-30s %7.4f %9.4f\n", MOMENT_NAMES[j], mm[j], dm[j]))
cat(sprintf("Mean inventory %.2f units; Pr(stock-out week) %.4f\n",
            mean(big$INV), mean(big$INV == 0 & big$D == 0)))

# 3. Three-step estimator on one panel
H <- 2000; T <- 104
panel <- simulate_panel(m, P0, TRUTH["alpha"], TRUTH["xi_B"], H, T, seed = 20261126)
s1 <- estimate_step1(m, panel)
s2 <- estimate_step2(m, panel, s1$theta)
counts <- build_counts(m, panel)
mle  <- estimate_step3_mle(m, counts, s2$omega, s2$Pi, spec = "common")
gmm  <- estimate_step3_gmm(m, panel, s2$omega, s2$Pi, spec = "common")
free <- estimate_step3_mle(m, counts, s2$omega, s2$Pi, spec = "free")

cat(sprintf("\nOne panel: H = %d, T = %d, %d purchases\n", H, T, s1$n_purchases))
cat("Step 1 (conditional logit, MLE)\n")
for (p in names(s1$theta)) cat(sprintf("  %-6s true %6.3f  est %6.3f (se %.3f)\n", p, TRUTH[p], s1$theta[p], s1$se[p]))
cat("Step 3 with one purchase-occasion utility (eta_small = eta_large)\n")
for (p in c("eta", "c")) cat(sprintf("  %-6s true %6.3f  MLE %6.3f (se %.3f)  GMM %6.3f (se %.3f)\n",
    p, TRUTH[p], mle$theta[p], mle$se[p], gmm$theta[p], gmm$se[p]))
cat(sprintf("  GMM Hansen J = %.2f (df %d, p = %.3f)\n", gmm$J, gmm$J_df, gmm$J_pvalue))
cat("Step 3 with size-specific utilities left free\n")
tr <- c(eta_small = 0.3, eta_large = 0.3, c = 0.04)
for (p in names(free$theta)) cat(sprintf("  %-10s true %6.3f  MLE %6.3f (se %.3f)\n", p, tr[p], free$theta[p], free$se[p]))
cat("  correlation of the free estimates:\n"); print(round(cov2cor(free$V), 3))
