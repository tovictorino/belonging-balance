# GTRE_ZH: persistent and transient inefficiency with determinants on both margins

R code for the four-component stochastic frontier used in *Belonging Has a Balance Sheet: System Affiliation, Output Definition, and Efficiency in Credit Cooperatives* (Victorino, submitted to the *Journal of Productivity Analysis*).

| File | Content |
|---|---|
| `psfm_gtre_zh.R` | The estimator: maximum simulated likelihood with Halton draws, analytic gradient, efficiency scores. |
| `psfm_gtre_zh_pf.R` | The same estimator with the option `draws = "perfirm"`, which gives each firm its own draws (as in sfa 1.2.0). |
| `demo_gtre_zh.R` | Simulated example. |
| `df_main.parquet` | Estimation panel: 6,902 observations of 864 Brazilian credit cooperatives, 2015–2024, values in 2024 reais. Built from public Central Bank of Brazil and Receita Federal data. |

Requires the R packages minqa, mvtnorm and randtoolbox, plus arrow to read the panel.

## Main model of the paper

```r
source("psfm_gtre_zh.R")
d <- as.data.frame(arrow::read_parquet("df_main.parquet"))
d <- transform(d, ly1 = -log(y1), lx1 = log(x1), lx2 = log(x2), lx3 = log(x3))
v <- c("lx1", "lx2", "lx3", "logAt", "li", "equity_ratio", "deposit_ratio", "age")
d[v] <- lapply(d[v], function(x) as.numeric(scale(x)))
for (y in 2016:2024) d[[paste0("year_", y)]] <- as.numeric(d$year_index == y)
f <- as.formula(paste("ly1 ~ lx1 + lx2 + lx3 +", paste0("year_", 2016:2024, collapse = " + "),
                      "| logAt + li + equity_ratio + filiacao + deposit_ratio | filiacao + age"))
fit <- psfm_gtre_zh(f, d, "CNPJ", inefdec = FALSE, halton_num = 1000, rand.gtre = 123,
                    maxit.bobyqa = 0, Method = "BFGS", maxit.optim = 3000)
```

The paper reports the highest likelihood (−654.12) over five starts: the default start and four perturbations of its optimum. Affiliation and age enter the persistent margin as cooperative means, hence the warning.

Code: GPL (≥ 2), derived from the R package sfa. Data: CC BY 4.0.
