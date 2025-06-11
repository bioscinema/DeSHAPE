# DADE  
**Distributional Analysis of Diversity Effects**

## Overview

**DADE** (*Distributional Analysis of Diversity Effects*) is an R package that provides a diagnostic framework for testing **distributional shifts** in alpha diversity across biological groups. It focuses on detecting shifts in:

- **Center** (e.g., median)
- **Dispersion** (spread)
- **Asymmetry** (skewness)

DADE supports both **permutation-based** and **quantile-regression-based** testing strategies, with or without covariate adjustment.

---

## Installation

Install the development version directly from GitHub:
```r
# install.packages("devtools")
devtools::install_github("bioscinema/DADE")
```

## Functions

| Function               | Description                                         | Best Used When                          |
|------------------------|-----------------------------------------------------|-----------------------------------------|
| `dade_perm_pair()`     | Permutation test for **2 groups**                   | Exploratory tests (center, spread, skew) |
| `dade_perm_multi()`    | Permutation test for **>2 groups**                  | Center or dispersion differences         |
| `wald_contrast_test()` | Quantile-regression-based **contrast test**         | Covariate-adjusted dispersion or skew    |


## 1. Permutation-Based Testing

### 1.1 `dade_perm_pair()`

Permutation-based comparison for **two groups** (e.g., group A vs group B).

```r
dade_perm_pair(response ~ group,
               data        = df,
               mode        = "center",     # or "dispersion", "skewness"
               alternative = "two.sided",  # only for "center" or "skewness"
               perm        = 999)
```

- `"center"`: Median difference test  
- `"dispersion"`: Median-based Levene test (Brown-Forsythe)  
- `"skewness"`: Quantile-asymmetry permutation test

---

### 1.2 `dade_perm_multi()`

Permutation-based comparison for **three or more groups**.

```r
dade_perm_multi(response ~ group,
                data = df,
                mode = "center",   # or "dispersion"
                perm = 999)
```

- `"center"`: Tests for median shift using permutation ANOVA  
- `"dispersion"`: Median-based Levene test

---

## 2. Quantile Regression Contrast Test

### 2.1 `wald_contrast_test()`

Wald-type test for linear contrasts of quantile regression coefficients across multiple quantile levels. Allows **confounder adjustment** and **tail-focused asymmetry testing**.

```r
wald_contrast_test(formula, data,
                   taus,           # e.g., c(0.25, 0.75)
                   contrast,       # numeric contrast vector
                   alternative = "two.sided",
                   kernel = "gaussian")
```

- `formula`: any valid R formula (e.g., `Shannon ~ group + Depth`)
- `taus`: quantiles of interest (e.g., `c(0.1, 0.5, 0.9)`)
- `contrast`: linear contrast vector applied to stacked coefficients
- `alternative`: `"two.sided"`, `"greater"`, or `"less"`

---

### 2.2 How to Construct the Contrast Vector

#### Step 1: Understand the model matrix layout

For a model `Shannon ~ group + covariate1 + covariate2`, with `K` quantiles, the full coefficient vector stacks like this:
```
(Intercept, group, cov1, cov2)[tau1],
(Intercept, group, cov1, cov2)[tau2],
...
(Intercept, group, cov1, cov2)[tauK]
```
Each quantile block contributes `(p + 1)` coefficients.

#### Step 2: Build a vector that contrasts coefficients across quantiles

- To test **dispersion**, subtract a covariate's effect at one quantile from its effect at another.
- To test **asymmetry**, apply weights like `+1, -2, +1` across the lower, median, and upper quantiles.

---

### 2.3 Examples (Illustrative Only – Do Not Copy)

#### Dispersion difference (adjusting nothing)

{r}
wald_contrast_test(Shannon ~ group_prefix,
                   data     = plot_df,
                   taus     = c(0.25, 0.75),
                   contrast = c(0, -1, 0, 1),
                   alternative = "greater")
{r\}

#### Dispersion difference (adjusting for sequencing depth)

{r}
wald_contrast_test(Shannon ~ group_prefix + Depth,
                   data     = plot_df,
                   taus     = c(0.25, 0.75),
                   contrast = c(0, -1, 0, 0, 1, 0),
                   alternative = "greater")
{r\}

#### Dispersion difference (adjusting for cohort and depth)

{r}
wald_contrast_test(Shannon ~ group_prefix + Cohort + Depth,
                   data     = plot_df,
                   taus     = c(0.25, 0.75),
                   contrast = c(0, -1, 0, 0, 0, 1, 0, 0),
                   alternative = "greater")
{r\}

#### Asymmetry difference (3 quantiles, unadjusted)

{r}
wald_contrast_test(Shannon ~ group_prefix,
                   data     = plot_df,
                   taus     = c(0.1, 0.5, 0.9),
                   contrast = c(0, 1, 0, -2, 0, 1),
                   alternative = "greater")
{r\}

#### Asymmetry difference (adjusting for cohort and depth)

{r}
wald_contrast_test(Shannon ~ group_prefix + Cohort + Depth,
                   data     = plot_df,
                   taus     = c(0.1, 0.5, 0.9),
                   contrast = c(0, 1, 0, 0, 0, -2, 0, 0, 0, 1, 0, 0),
                   alternative = "greater")
{r\}

---

## Workflow Summary

1. Use `dade_perm_pair()` or `dade_perm_multi()` for unadjusted, distribution-based testing.
2. Use `wald_contrast_test()` when:
   - You want to adjust for covariates (e.g., sequencing depth).
   - You want to assess tail asymmetry more precisely.
3. Build your `contrast` vector carefully by indexing the relevant coefficient positions across quantile blocks.

---
