# Monte Carlo results

True values: alpha = 2.00, xi_B = -0.30, eta = 0.30, c = 0.04.

## E1 Baseline, H = 500 households, T = 104 weeks (200 replications)

| estimator / parameter | mean | bias | RMSE | mean SE | 95% coverage |
|---|---|---|---|---|---|
| Step 1 MLE: alpha | 2.0054 | +0.0054 | 0.0379 | 0.0367 | 0.960 |
| Step 1 MLE: xi_B | -0.2993 | +0.0007 | 0.0157 | 0.0154 | 0.965 |
| Step 3 MLE: eta | 0.3013 | +0.0013 | 0.0314 | 0.0225 | 0.845 |
| Step 3 MLE: c | 0.0398 | -0.0002 | 0.0015 | 0.0012 | 0.880 |
| Step 3 MLE, true steps 1-2: eta | 0.2989 | -0.0011 | 0.0234 | 0.0225 | 0.940 |
| Step 3 MLE, true steps 1-2: c | 0.0399 | -0.0001 | 0.0012 | 0.0012 | 0.970 |
| Step 3 GMM: eta | 0.3019 | +0.0019 | 0.0334 | 0.0269 | 0.890 |
| Step 3 GMM: c | 0.0398 | -0.0002 | 0.0030 | 0.0030 | 0.950 |

Hansen J test rejection rate at 5%: 0.075.

## E2 Size-specific purchase utilities left free, H = 500

| parameter | true | mean | sd | RMSE |
|---|---|---|---|---|
| eta_small | 0.300 | 0.2751 | 0.4889 | 0.4883 |
| eta_large | 0.300 | 0.2488 | 0.9745 | 0.9734 |
| c | 0.040 | 0.0392 | 0.0120 | 0.0120 |
| eta_large - 2 eta_small | -0.300 | -0.3013 | 0.0314 | 0.0314 |

RMSE of c: 0.0015 with the restriction, 0.0120 without (ratio 7.9).

## E1 Baseline, H = 2000 households, T = 104 weeks (200 replications)

| estimator / parameter | mean | bias | RMSE | mean SE | 95% coverage |
|---|---|---|---|---|---|
| Step 1 MLE: alpha | 1.9995 | -0.0005 | 0.0186 | 0.0183 | 0.960 |
| Step 1 MLE: xi_B | -0.2996 | +0.0004 | 0.0079 | 0.0077 | 0.945 |
| Step 3 MLE: eta | 0.3007 | +0.0007 | 0.0152 | 0.0112 | 0.845 |
| Step 3 MLE: c | 0.0400 | -0.0000 | 0.0008 | 0.0006 | 0.895 |
| Step 3 MLE, true steps 1-2: eta | 0.3011 | +0.0011 | 0.0107 | 0.0112 | 0.945 |
| Step 3 MLE, true steps 1-2: c | 0.0400 | -0.0000 | 0.0006 | 0.0006 | 0.955 |
| Step 3 GMM: eta | 0.3003 | +0.0003 | 0.0149 | 0.0135 | 0.930 |
| Step 3 GMM: c | 0.0401 | +0.0001 | 0.0014 | 0.0015 | 0.965 |

Hansen J test rejection rate at 5%: 0.075.

## E2 Size-specific purchase utilities left free, H = 2000

| parameter | true | mean | sd | RMSE |
|---|---|---|---|---|
| eta_small | 0.300 | 0.2946 | 0.2366 | 0.2360 |
| eta_large | 0.300 | 0.2885 | 0.4708 | 0.4698 |
| c | 0.040 | 0.0398 | 0.0059 | 0.0059 |
| eta_large - 2 eta_small | -0.300 | -0.3007 | 0.0152 | 0.0152 |

RMSE of c: 0.0008 with the restriction, 0.0059 without (ratio 7.7).

## E4 Panel composition and the initial-inventory prior (N x T = 208,000)

| H x T | prior | bias eta | RMSE eta | bias c | RMSE c | coverage eta | coverage c |
|---|---|---|---|---|---|---|---|
| 8000 x 26 | stationary (correct) | -0.0033 | 0.0148 | +0.0001 | 0.0008 | 0.890 | 0.890 |
| 8000 x 26 | empty (wrong) | -1.1160 | 1.1161 | +0.0628 | 0.0628 | 0.000 | 0.000 |
| 4000 x 52 | stationary (correct) | -0.0001 | 0.0155 | +0.0001 | 0.0008 | 0.800 | 0.890 |
| 4000 x 52 | empty (wrong) | -0.5627 | 0.5629 | +0.0241 | 0.0241 | 0.000 | 0.000 |
| 2000 x 104 | stationary (correct) | +0.0022 | 0.0147 | -0.0001 | 0.0008 | 0.870 | 0.900 |
| 2000 x 104 | empty (wrong) | -0.2787 | 0.2791 | +0.0101 | 0.0101 | 0.000 | 0.000 |
| 1000 x 208 | stationary (correct) | +0.0007 | 0.0159 | -0.0000 | 0.0008 | 0.840 | 0.910 |
| 1000 x 208 | empty (wrong) | -0.1399 | 0.1409 | +0.0046 | 0.0047 | 0.000 | 0.000 |
