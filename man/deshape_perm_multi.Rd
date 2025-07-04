% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dashape_perm_multi.R
\name{deshape_perm_multi}
\alias{deshape_perm_multi}
\title{Permutation-Based Multi-Group Comparison}
\usage{
deshape_perm_multi(
  formula,
  data,
  mode = c("center", "dispersion", "skewness"),
  perm = 999,
  seed
)
}
\arguments{
\item{formula}{A formula of the form `response ~ group`, where the group variable has more than two levels.}

\item{data}{A data frame containing the variables in the formula.}

\item{mode}{Type of test to perform. Must be one of `"center"` (permutation median ANOVA), `"dispersion"` (IQR-based permutation test), or `"skewness"` (quantile-based asymmetry test).}

\item{perm}{Number of permutations. Default is 999. Used for all modes.}

\item{seed}{Random seed for reproducibility. Default is NULL.}
}
\value{
This function prints the permutation p-value to the console.
}
\description{
Performs a permutation-based test across multiple groups to assess differences in
center (median), dispersion (IQR-based), or asymmetry (quantile-based).
}
\examples{
\dontrun{
data <- data.frame(Shannon = rnorm(90), Group = rep(c("A", "B", "C"), each = 30))
deshape_perm_multi(Shannon ~ Group, data = data, mode = "center")
}
}
