## Simulated example: estimate GTRE_ZH, compare with the true parameters and
## efficiencies, and test the persistent determinant by likelihood ratio.

source("psfm_gtre_zh.R")
set.seed(42)

N <- 120; T <- 6
beta_true  <- c(cons = 0.5, x1 = 0.5, x2 = 0.5)
sigv_true  <- 0.15
sigr_true  <- 0.20
delta_true <- c(d0 = log(0.04), d1 = 0.8)
gamma_true <- c(g0 = log(0.09), g1 = 0.9)

id   <- rep(1:N, each = T)
x1   <- rnorm(N * T)
x2   <- rnorm(N * T)
zu   <- rnorm(N * T)
zh_i <- rnorm(N)
zh   <- zh_i[id]

sig_u <- sqrt(exp(delta_true[1] + delta_true[2] * zu))
sig_h <- sqrt(exp(gamma_true[1] + gamma_true[2] * zh_i))
r_i   <- rnorm(N, 0, sigr_true)
eta_i <- abs(rnorm(N, 0, 1)) * sig_h
v_it  <- rnorm(N * T, 0, sigv_true)
u_it  <- abs(rnorm(N * T)) * sig_u

y <- beta_true[1] + beta_true[2] * x1 + beta_true[3] * x2 +
     r_i[id] - eta_i[id] + v_it - u_it
dat <- data.frame(name = id, y = y, x1 = x1, x2 = x2, zu = zu, zh = zh)

fit <- psfm_gtre_zh(y ~ x1 + x2 | zu | zh, data = dat, individual = "name",
                    halton_num = 150, rand.gtre = 123, compute_te = TRUE)
print(fit)

truth <- c(sigv_true, sigr_true, beta_true, delta_true, gamma_true)
print(round(cbind(true = unname(truth), estimate = fit$coefficients), 3))

print(round(c(transient  = cor(fit$U, exp(-u_it)),
              persistent = cor(fit$H, exp(-eta_i)),
              overall    = cor(fit$O, exp(-(eta_i[id] + u_it)))), 3))

fit0 <- psfm_gtre_zh(y ~ x1 + x2 | zu | 1, data = dat, individual = "name",
                     halton_num = 150, rand.gtre = 123, maxit.bobyqa = 2000,
                     maxit.optim = 500, compute_te = FALSE)
LR <- 2 * (fit$loglik - fit0$loglik)
df <- length(fit$coefficients) - length(fit0$coefficients)
cat("LR =", round(LR, 2), " df =", df, " p =",
    format.pval(pchisq(LR, df, lower.tail = FALSE)), "\n")

print(lapply(margeff_gtre_zh(fit)[1:2], round, 4))
