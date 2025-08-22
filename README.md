# DeSHAPE  
**Decomposing ecological Structure through Heterogeneity, Asymmetry and Pattern Evaluation of alpha-Diversity**

## Overview

**DeSHAPE** (*Decomposing ecological Structure through Heterogeneity, Asymmetry and Pattern Evaluation of alpha-Diversity*) is an R package that provides a diagnostic framework for testing **distributional shifts** in alpha diversity across biological groups. It focuses on detecting shifts in:

- **Center** (e.g., median)
- **Dispersion** (spread and variability)
- **Asymmetry** (skewness)

DeSHAPE supports both **permutation-based** and **quantile-regression-based** testing strategies, with or without covariate adjustment.

---

## Installation

Install the development version directly from GitHub:
```r
# install.packages("devtools")
library(devtools)
devtools::install_github("bioscinema/DeSHAPE")
library(DeSHAPE)
```

## Functions

| Function                  | Description                                               | Best Used When                                                             |
|---------------------------|-----------------------------------------------------------|----------------------------------------------------------------------------|
| `deshape_perm_pair()`     | Permutation test for **2 groups**                         | Exploratory tests (center, dispersion, asymmetry) for two groups                                 |
| `deshape_perm_multi()`    | Permutation test for **>2 groups**                        | Exploratory tests (center, dispersion, asymmetry) for multiple groups                                                |
| `deshape_wald_contrast()` | Quantile-regression-based **contrast test**               | Covariate-adjusted tests for center, dispersion, or asymmetry (sample size of each group > 30)                            |
| `deshape_glm_resid_test()`| Residual-shape permutation test **after** a GLM fit       | Covariate-adjusted tests for center, dispersion, or asymmetry   |




## 1. Permutation-Based Testing

### 1.1 `deshape_perm_pair()`

Permutation-based comparison for **two groups** (e.g., group A vs group B).

```r
deshape_perm_pair(response ~ group,
               data        = df,
               mode        = "center",     # or "dispersion", "skewness"
               alternative = "two.sided",  
               perm        = 999,
               seed        = 2025)
```

- `"center"`: Median difference permutation test  
- `"dispersion"`: Interquartile range (IQR) permutation test
- `"skewness"`: Quantile-based asymmetry score permutation test

---

### 1.2 `deshape_perm_multi()`

Permutation-based comparison for **three or more groups** (center, dispersion, or skewness).

```r
deshape_perm_multi(response ~ group,
                data = df,
                mode = "center",   # or "dispersion", "skewness"
                perm = 999,
                seed = 2025)
```

- `"center"`: Tests for median shift using permutation test on ANOVA-style test statistic  
- `"dispersion"`: Tests for IQR shift using permutation test on ANOVA-style test statistic  
- `"skewness"`: Tests for asymmetry score shift using permutation test on ANOVA-style test statistic  

---

\textbf{2. Quantile Regression Contrast Test}

\textbf{2.1 \texttt{deshape\_wald\_contrast()}}

\texttt{deshape\_wald\_contrast()} performs a Wald-type test on linear contrasts 
of quantile-regression coefficients across multiple quantile levels. 
It allows both \textbf{pre-specified contrasts} 
(\texttt{mode = "customize"}) and \textbf{automatically constructed contrasts} 
for common hypotheses about \textbf{dispersion} and \textbf{symmetry}.

The function fits separate quantile regressions at each specified quantile 
(\texttt{taus}) using \texttt{quantreg::rq}, stacks the coefficients, 
and builds an asymptotic covariance matrix using kernel density estimates 
of residuals at 0. A Wald $\chi^2$ statistic (1 d.f.) is then computed, 
and p-values are returned under the chosen alternative hypothesis.

\begin{verbatim}
deshape_wald_contrast(
  formula, data,
  mode = c("customize", "dispersion", "symmetry"),
  taus = NULL,
  contrast = NULL,
  alternative = c("two.sided", "greater", "less"),
  kernel = "gaussian"
)
\end{verbatim}

\textbf{Arguments}
\begin{itemize}
  \item \texttt{formula}: model formula of the form \texttt{response ~ predictors}. 
        For non-custom modes (\texttt{"dispersion"}, \texttt{"symmetry"}), 
        the first non-intercept term must be a binary group variable 
        that produces exactly one column in \texttt{model.matrix}.
  \item \texttt{data}: a data frame containing all variables in the formula.
  \item \texttt{mode}: how to construct the contrast.
    \begin{itemize}
      \item \texttt{"dispersion"} – automatically compares group effects 
            between two quantiles (default \texttt{taus = c(0.25, 0.75)}).
      \item \texttt{"symmetry"} – automatically compares group effects 
            across three quantiles with weights +1, –2, +1 
            (default \texttt{taus = c(0.1, 0.5, 0.9)}).
      \item \texttt{"customize"} – user must provide both \texttt{taus} 
            and a full contrast vector.
    \end{itemize}
  \item \texttt{taus}: numeric vector of quantile levels in ascending order. 
        Required for \texttt{mode = "customize"}. Must have length 2 
        for \texttt{"dispersion"} and 3 for \texttt{"symmetry"}.
  \item \texttt{contrast}: numeric vector (or single-column matrix) of 
        length $p \times K$, where $p = \texttt{ncol(model.matrix(...))}$ 
        and $K = \texttt{length(taus)}$. Only used when 
        \texttt{mode = "customize"}.
  \item \texttt{alternative}: \texttt{"two.sided"} (default), 
        \texttt{"greater"}, or \texttt{"less"}.
  \item \texttt{kernel}: currently only \texttt{"gaussian"} is supported.
\end{itemize}

\textbf{Value}
\begin{itemize}
  \item \texttt{test\_stat}: Wald $\chi^2$ statistic (df = 1).
  \item \texttt{p\_value}: p-value computed using the chosen alternative hypothesis.
\end{itemize}

\textbf{2.2 Coefficient stacking and contrast construction}

At each quantile $\tau$, the regression has $p$ coefficients. 
These are \textbf{stacked by quantile} into a vector of length $p \times K$:

\[
[\beta_{\text{Int}}(\tau_1), \beta_{\text{group}}(\tau_1), 
 \beta_{\text{cov1}}(\tau_1), \ldots, 
 \beta_{\text{Int}}(\tau_K), \beta_{\text{group}}(\tau_K), 
 \beta_{\text{cov1}}(\tau_K), \ldots ]
\]

A \textbf{contrast vector} selects and compares coefficients across quantiles.
\begin{itemize}
  \item \textbf{Dispersion test (auto):} with $\tau = (\tau_L, \tau_U)$, 
        places –1 on the group coefficient at $\tau_L$ 
        and +1 at $\tau_U$.
  \item \textbf{Symmetry test (auto):} with 
        $\tau = (\tau_L, \tau_M, \tau_U)$, 
        places +1, –2, +1 on the group coefficients 
        at $\tau_L$, $\tau_M$, and $\tau_U$.
  \item \textbf{Customize:} any valid contrast can be supplied 
        (length = $p \times K$).
\end{itemize}

\textbf{2.3 Examples}

\textbf{Example A: Dispersion test (automatic)}  
Tests whether the group effect increases from $\tau = 0.25$ to $\tau = 0.75$.
\begin{verbatim}
deshape_wald_contrast(
  Shannon ~ group + Depth,
  data = df,
  mode = "dispersion"
)
\end{verbatim}

\textbf{Example B: Symmetry test (automatic)}  
Tests whether the group effect is symmetric across $\tau = 0.1, 0.5, 0.9$.
\begin{verbatim}
deshape_wald_contrast(
  Shannon ~ group + Depth,
  data = df,
  mode = "symmetry"
)
\end{verbatim}

\textbf{Example C: Custom contrast}  
Explicitly test whether $\beta_{\text{group}}(0.75) - \beta_{\text{group}}(0.25) = 0$. 
Here, $p = 3$ (Intercept, group, Depth) and $K = 2$ 
($\texttt{taus} = 0.25, 0.75$), so $\texttt{length(contrast)} = 6$.
\begin{verbatim}
taus <- c(0.25, 0.75)
contrast <- c(0, -1, 0, 0, +1, 0)

deshape_wald_contrast(
  Shannon ~ group + Depth,
  data = df,
  mode = "customize",
  taus = taus,
  contrast = contrast,
  alternative = "greater"
)
\end{verbatim}

\textbf{2.4 Practical notes}
\begin{itemize}
  \item For median differences only, a single quantile regression 
        at $\tau = 0.5$ followed by \texttt{summary()} 
        may be simpler and more efficient.
  \item If the group variable is not binary or expands to multiple columns 
        in \texttt{model.matrix}, you must use 
        \texttt{mode = "customize"}.
  \item Ensure \texttt{taus} are strictly ascending.
\end{itemize}


## 3. GLM Residual Shape Test

### 3.1 `deshape_glm_resid_test()`

Compares the **center**, **dispersion**, or **skewness** of GLM residuals
across user-defined groups while automatically adjusting for all other
covariates in the model.

```r
# adjust for Depth, Cohort, then compare residual dispersion by group_prefix
deshape_glm_resid_test(Shannon ~ group_prefix + Cohort + Depth,
                  data = plot_df,
                  mode = "center",
                  family = Gamma(link = "log"), 
                  group = "group_prefix", 
                  alternative = "greater",
                  q = 0.25, 
                  B = 999, 
                  seed = 2025)
```

Internally the function drops `group` from the right-hand side before
fitting the GLM, so residuals are free of the group effect yet still
reflect all other covariate adjustments.

**Choosing `q`**

- **`mode = "center"`** – `q` is ignored; you may supply any value.
- **`mode = "dispersion"`** – set **`q = 0.25`** to compute the classical inter-quartile range `Q(0.75) - Q(0.25)`.
- **`mode = "skewness"`** – `q` controls the tail width in the asymmetry score:  
  - **`q = 0.25`** → `Q(0.75) - 2 × Q(0.5) + Q(0.25)`  
  - **`q = 0.10`** → `Q(0.90) - 2 × Q(0.5) + Q(0.10)`


---

## Workflow Summary

1. Use `deshape_perm_pair()` or `deshape_perm_multi()` for unadjusted, distribution-based testing.
2. Use `deshape_wald_contrast()` when you want formal quantile-regression
   contrasts (especially for dispersion or tail asymmetry) **and** need
   covariate adjustment; for a simple median shift a single τ = 0.5
   quantile regression followed by `summary()` is usually quicker. Build your `contrast` vector carefully by indexing the relevant coefficient positions across quantile blocks.
3. Use `deshape_glm_resid_test()` when you have limited sample size or want it as sensitivity/robustness check to Wald contrast method.

---
