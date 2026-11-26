# True parameter values used to generate all simulated data
TRUTH <- c(alpha = 2.0,   # price sensitivity
           xi_B = -0.3,   # preference for brand B relative to brand A
           eta = 0.3,     # utility of a purchase occasion (same for both sizes)
           c = 0.04)      # storage cost per unit per week
TRUE_THETA3 <- c(TRUTH["eta"], TRUTH["eta"], TRUTH["c"])
