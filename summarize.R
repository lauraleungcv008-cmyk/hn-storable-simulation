# Tables and figures from the Monte Carlo output in results/
#   Rscript summarize.R
source("R/truth.R")
stat <- function(est, se, truth) {
  c(mean = mean(est), bias = mean(est) - truth, rmse = sqrt(mean((est - truth)^2)),
    mean_se = mean(se), coverage = mean(abs(est - truth) <= 1.96 * se, na.rm = TRUE))
}
row <- function(label, s) sprintf("| %s | %.4f | %+.4f | %.4f | %.4f | %.3f |", label, s[1], s[2], s[3], s[4], s[5])
head5 <- c("| estimator / parameter | mean | bias | RMSE | mean SE | 95% coverage |", "|---|---|---|---|---|---|")
out <- c("# Monte Carlo results", "",
         sprintf("True values: alpha = %.2f, xi_B = %.2f, eta = %.2f, c = %.2f.", TRUTH[1], TRUTH[2], TRUTH[3], TRUTH[4]))

# E1 and E2 ---------------------------------------------------------------
for (H in c(500, 2000)) {
  r <- readRDS(sprintf("results/baseline_H%d_T104.rds", H))
  out <- c(out, "", sprintf("## E1 Baseline, H = %d households, T = 104 weeks (%d replications)", H, nrow(r)), "", head5,
    row("Step 1 MLE: alpha", stat(r[, "s1.alpha"], r[, "s1se.alpha"], TRUTH["alpha"])),
    row("Step 1 MLE: xi_B", stat(r[, "s1.xi_B"], r[, "s1se.xi_B"], TRUTH["xi_B"])),
    row("Step 3 MLE: eta", stat(r[, "mle.eta"], r[, "mlese.eta"], TRUTH["eta"])),
    row("Step 3 MLE: c", stat(r[, "mle.c"], r[, "mlese.c"], TRUTH["c"])),
    row("Step 3 MLE, true steps 1-2: eta", stat(r[, "orc.eta"], r[, "orcse.eta"], TRUTH["eta"])),
    row("Step 3 MLE, true steps 1-2: c", stat(r[, "orc.c"], r[, "orcse.c"], TRUTH["c"])),
    row("Step 3 GMM: eta", stat(r[, "gmm.eta"], r[, "gmmse.eta"], TRUTH["eta"])),
    row("Step 3 GMM: c", stat(r[, "gmm.c"], r[, "gmmse.c"], TRUTH["c"])),
    "", sprintf("Hansen J test rejection rate at 5%%: %.3f.", mean(r[, "Jp"] < 0.05)),
    "", sprintf("## E2 Size-specific purchase utilities left free, H = %d", H), "",
    "| parameter | true | mean | sd | RMSE |", "|---|---|---|---|---|")
  fr <- list(c("free.eta_small", 0.3), c("free.eta_large", 0.3), c("free.c", 0.04), c("contrast", -0.3))
  lab <- c("eta_small", "eta_large", "c", "eta_large - 2 eta_small")
  for (j in 1:4) { x <- r[, fr[[j]][1]]; tv <- as.numeric(fr[[j]][2])
    out <- c(out, sprintf("| %s | %.3f | %.4f | %.4f | %.4f |", lab[j], tv, mean(x), sd(x), sqrt(mean((x - tv)^2)))) }
  out <- c(out, "", sprintf("RMSE of c: %.4f with the restriction, %.4f without (ratio %.1f).",
      sqrt(mean((r[, "mle.c"] - 0.04)^2)), sqrt(mean((r[, "free.c"] - 0.04)^2)),
      sqrt(mean((r[, "free.c"] - 0.04)^2)) / sqrt(mean((r[, "mle.c"] - 0.04)^2))))
}

# E4 ----------------------------------------------------------------------
designs <- list(c(8000, 26), c(4000, 52), c(2000, 104), c(1000, 208))
out <- c(out, "", "## E4 Panel composition and the initial-inventory prior (N x T = 208,000)", "",
         "| H x T | prior | bias eta | RMSE eta | bias c | RMSE c | coverage eta | coverage c |", "|---|---|---|---|---|---|---|---|")
E4 <- NULL
for (d in designs) {
  f <- sprintf("results/composition_H%d_T%d.rds", d[1], d[2]); if (!file.exists(f)) next
  r <- readRDS(f)
  for (p in c("stat", "empty")) {
    e <- stat(r[, paste0(p, ".eta")], r[, paste0(p, "se.eta")], TRUTH["eta"])
    cc <- stat(r[, paste0(p, ".c")], r[, paste0(p, "se.c")], TRUTH["c"])
    out <- c(out, sprintf("| %d x %d | %s | %+.4f | %.4f | %+.4f | %.4f | %.3f | %.3f |", d[1], d[2],
        ifelse(p == "stat", "stationary (correct)", "empty (wrong)"), e[2], e[3], cc[2], cc[3], e[5], cc[5]))
    E4 <- rbind(E4, data.frame(T = d[2], prior = p, bias_eta = e[2], rmse_eta = e[3], bias_c = cc[2], rmse_c = cc[3]))
  }
}
writeLines(out, "results/montecarlo_tables.md"); cat(out, sep = "\n")

# Figures -----------------------------------------------------------------
r <- readRDS("results/baseline_H2000_T104.rds")
png("figures/fig1_identification.png", 2400, 900, res = 220)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1), cex = 0.8)
plot(r[, "free.eta_small"], r[, "free.c"], pch = 16, col = "#4C6A9280", xlab = "estimate of eta_small",
     ylab = "estimate of c", main = "Free size utilities: estimates lie on a ridge")
points(0.3, 0.04, pch = 4, col = "#C8553D", cex = 2, lwd = 3)
plot(r[, "mle.eta"], r[, "mle.c"], pch = 16, col = "#4C6A9280", xlab = "estimate of eta",
     ylab = "estimate of c", main = "With eta_small = eta_large")
points(0.3, 0.04, pch = 4, col = "#C8553D", cex = 2, lwd = 3)
legend("topright", legend = "true value", pch = 4, col = "#C8553D", pt.lwd = 3, bty = "n")
dev.off()

if (!is.null(E4)) {
  png("figures/fig2_composition.png", 2400, 900, res = 220)
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1), cex = 0.8)
  for (v in c("eta", "c")) {
    b <- E4[[paste0("bias_", v)]] / TRUTH[v] * 100
    plot(NA, xlim = range(log2(E4$T)), ylim = range(c(b, 0)), xaxt = "n", xlab = "weeks per household (N x T fixed)",
         ylab = "bias, % of true value", main = paste("Bias of", v))
    axis(1, at = log2(unique(E4$T)), labels = unique(E4$T)); abline(h = 0, lty = 3)
    lines(log2(E4$T[E4$prior == "stat"]), b[E4$prior == "stat"], type = "b", pch = 16, col = "#4C6A92", lwd = 2)
    lines(log2(E4$T[E4$prior == "empty"]), b[E4$prior == "empty"], type = "b", pch = 15, col = "#C8553D", lwd = 2)
    if (v == "eta") legend("bottomright", c("correct initial-inventory prior", "wrong prior (all start empty)"),
                           col = c("#4C6A92", "#C8553D"), pch = c(16, 15), lwd = 2, bty = "n")
  }
  dev.off()
}
