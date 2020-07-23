Example PMI Calculation
================
Zachary del Rosario
2020-07-09

This companion notebook accompanies the manuscript “Precision Materials
Indices: Materials Selection with Statistically-rigorous Reliability
Analysis.” It provides a step-by-step example of computing a precision
materials index (PMI) from raw data.

# Setup

<!-- -------------------------------------------------- -->

The following is just setup for the
    example.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✔ ggplot2 3.3.1     ✔ purrr   0.3.4
    ## ✔ tibble  3.0.1     ✔ dplyr   1.0.0
    ## ✔ tidyr   1.1.0     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
## Ground-truth parameters for Al 6061-T6
median_Y <- 40000.0
median_rho <- 0.0975
cov_Y <- 0.03
cov_rho <- 0.01

## ## Ground-truth parameters for Ti 6Al 4V
## median_Y <- 160000.0
## median_rho <- 0.1600
## cov_Y <- 0.06
## cov_rho <- 0.02

## Derived parameters
mu_Y <- log(median_Y)
mu_rho <- log(median_rho)
sig_Y <- cov_Y
sig_rho <- cov_rho

meanlog_Y <- mu_Y + 0.5 * sig_Y^2
meanlog_rho <- mu_rho + 0.5 * sig_rho^2
sdlog_Y <- cov_Y
sdlog_rho <- cov_rho
```

# Example Calculation

<!-- -------------------------------------------------- -->

The following are the steps of calculating a PMI.

## 1\. Estimate parameters

<!-- ------------------------- -->

For demonstration purposes, simulate some data.

``` r
## Sample size
n <- 10

## Generate synthetic data
set.seed(101) # For reproducibility
df_data <-
  tibble(
    Y = rlnorm(n, meanlog = meanlog_Y, sdlog = sdlog_Y),
    rho = rlnorm(n, meanlog = meanlog_rho, sdlog = sdlog_rho)
  )
```

Estimate parameters from simulated data.

``` r
## KEY STEP: Log-transform data
df_transformed <-
  df_data %>%
  mutate(log_rho = log(rho), log_Y = log(Y))

## Estimate parameters
mu_hat <-
  df_transformed %>%
  summarize(across(c(log_rho, log_Y), mean)) %>%
  unlist(.)

sigma_hat_11 <- var(df_transformed$log_rho)
sigma_hat_22 <- var(df_transformed$log_Y)
sigma_hat_12 <- cov(
  x = df_transformed$log_rho,
  y = df_transformed$log_Y
)

Sigma_hat <- matrix(
  c(sigma_hat_11, sigma_hat_12,
    sigma_hat_12, sigma_hat_22),
  nrow = 2
)
```

## 2\. Convert to monomial-standard form

<!-- ------------------------- -->

In monomial-standard form the beam strength cost is

\[\min_{t}\, Q = (l) t^{p=2} \rho^{v_{\rho} = 1},\]

while the reliability constraint
is

\[\text{s.t. } \mathbb{P}_{Y}[ t^{q=-3} Y^{w_Y=-1} \leq (6fl)^{-1}] \geq \mathcal{R}.\]

## 3\. Extract and compute

<!-- ------------------------- -->

The extracted quantities are

``` r
p <- 2
q <- -3
v <- c("rho" = 1, "Y" = 0)
w <- c("rho" = 0, "Y" = -1)
```

therefore the estimated statistics are

``` r
mu_v <- v %*% mu_hat
mu_w <- w %*% mu_hat
sig_w <- sqrt(t(w) %*% Sigma_hat %*% w)
```

## 4\. Set target and compute MIB

<!-- ------------------------- -->

``` r
## Set targets
beta_t <- 3 # Reliability index
C <- 0.95 # Confidence level

## Compute the MIB
beta <- beta_t # Aspirational approximation
b <- qt(p = C, ncp = sqrt(n) * beta, df = n) / sqrt(n) - beta
```

## 5\. Compute the PMI factors

<!-- ------------------------- -->

``` r
MMI <- exp(p/q * mu_w - mu_v)
f_beta_t <- exp(p/q * beta_t * sig_w)
f_b <- exp(p/q * b * sig_w)

PMI <- MMI * f_beta_t * f_b
knockdown <- 1 - f_beta_t * f_b

tibble(MMI = MMI, PMI = PMI, f_beta_t = f_beta_t, f_b = f_b) %>%
  knitr::kable(digits = 2)
```

|      MMI |      PMI | f\_beta\_t |     f\_b |
| -------: | -------: | ---------: | -------: |
| 12111.05 | 11435.99 |  0.9655435 | 0.977958 |

``` r
knockdown
```

    ##            [,1]
    ## [1,] 0.05573899

Variability (COV) of the material strength is modest: The combined
knockdown of the material index due to aleatory and epistemic factors is
about \(6\%\).
