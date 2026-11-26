# Hendel-Nevo (2006) three-step estimator: replication on a small model and a simulation study

Author: Shihui Liang. Language: R (base R only, plus `parallel`).

## What it does
1. **Model** (`R/model.R`, `R/simulate.R`): a storable-good model with 2 brands x 2 pack sizes,
   Markov sale prices, unit weekly consumption and unobserved household inventory.
   Solved by value function iteration and by policy iteration (they agree to 1e-11).
2. **Estimator** (`R/step1_*.R`, `R/step2_*.R`, `R/step3_mle.R`): the three steps of Hendel and Nevo (2006).
   Step 1: brand choice conditional on size, conditional logit by MLE.
   Step 2: inclusive values per size and their Markov transition.
   Step 3: dynamic size choice by nested fixed point MLE, integrating over unobserved initial inventory.
   `R/step3_gmm.R` re-estimates step 3 by efficient GMM on six moments.
3. **Simulation study** (`run_montecarlo.R`, `summarize.R`):
   E1 baseline validation; E2 a restriction on purchase utility; E4 panel composition
   and the initial-inventory assumption.

## Reproduce
```
Rscript run_single.R          # checks and one estimation, ~10 seconds
Rscript run_montecarlo.R      # full study, ~45 minutes on 2 cores
Rscript summarize.R           # results/montecarlo_tables.md and figures/
```
Seeds are fixed, so results reproduce exactly.

## Simplifications relative to Hendel and Nevo (2006)
Two brands, discrete price states, fixed consumption of one unit per week, discrete inventory,
homogeneous households. The discount factor (0.95) and the stock-out cost are fixed, not estimated.
Step 3 standard errors do not include the Murphy-Topel correction for steps 1-2.
