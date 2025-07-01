#' DeSHAPE GLM‐Residual Group Shape Test
#'
#' Implements the **DeSHAPE** framework on residuals from a
#' generalized-linear model (GLM) to compare *shape* characteristics
#' (location, dispersion, or skewness) of the residual distribution
#' across user-defined groups.  The procedure
#' \enumerate{
#'   \item fits a GLM given by \code{formula} and extracts Pearson residuals;
#'   \item computes a group-specific summary statistic that targets the
#'         chosen aspect of shape (\code{"center"}, \code{"dispersion"},
#'         or \code{"skewness"});
#'   \item forms a weighted cross-group contrast statistic, \eqn{T};
#'   \item obtains a permutation reference distribution of \eqn{T} by
#'         randomly re-shuffling group labels \code{B} times;
#'   \item returns the observed statistic and a permutation
#'         \eqn{p}-value.
#' }
#' The test generalises the original DeSHAPE methodology to any GLM
#' supported by \code{\link[stats]{glm}} (e.g., Gaussian, Poisson,
#' binomial) and is especially useful when covariates must be adjusted
#' before shape differences are interrogated.
#'
#' @section Test statistics:
#' \describe{
#'   \item{\code{mode = "center"}}{Median of residuals in each group.}
#'   \item{\code{mode = "dispersion"}}{Inter-quantile range
#'         \eqn{Q_{1-q} - Q_q}.}
#'   \item{\code{mode = "skewness"}}{Third-order, quantile-based
#'         contrast \eqn{Q_{1-q} - 2Q_{0.5} + Q_q}.}
#' }
#'
#' @param formula A model \code{\link[stats]{formula}} passed to
#'   \code{\link[stats]{glm}}.
#' @param data    A data frame containing the variables referenced in
#'   \code{formula} and the grouping variable.
#' @param mode    Character; which aspect of residual shape to test.
#'   Must be one of \code{"center"}, \code{"dispersion"},
#'   \code{"skewness"}.  See *Test statistics*.
#' @param family  A GLM family object; default is
#'   \code{\link[stats]{gaussian}}.
#' @param group   **Character string** giving the column name in
#'   \code{data} that defines groups whose residual shapes are
#'   compared.
#' @param alternative Alternative hypothesis: \code{"two.sided"}
#'   (default), \code{"greater"}, or \code{"less"}.
#' @param q       Lower tail probability used in the inter-quantile
#'   contrasts; must satisfy \eqn{0 < q < 0.5}.  Default \code{0.25}.
#' @param B       Integer; number of label permutations (default
#'   \code{999}).
#' @param seed    Optional integer seed for reproducibility.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{statistic}}{Observed test statistic \eqn{T}.}
#'   \item{\code{p.value}}{Permutation \eqn{p}-value.}
#'   \item{\code{T_perm}}{Vector of permuted statistics (length
#'         \code{B}).}
#'   \item{\code{mode}}{Requested shape aspect.}
#'   \item{\code{method}}{Character describing the test.}
#'   \item{\code{B}}{Number of permutations performed.}
#' }
#'
#' @details
#' The statistic aggregates group-level summaries through a
#' weighted mean‐of-squares so that groups with more observations
#' contribute proportionally.  Under the null hypothesis of equal
#' residual shape across groups, permuting group labels yields an
#' exact test when residuals are exchangeable and an
#' asymptotically valid test otherwise.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n  <- 120
#' g  <- gl(3, n/3, labels = c("A", "B", "C"))
#' x1 <- rnorm(n)
#' y  <- 0.5 + 1.2 * x1 + rnorm(n, sd = 0.8)          # Gaussian GLM
#' dat <- data.frame(y, x1, g)
#'
#' # Test equality of residual dispersion across groups
#' out <- deshape_glm_resid_test(y ~ x1,
#'                               data  = dat,
#'                               family = gaussian(),
#'                               group  = "g",
#'                               mode   = "dispersion",
#'                               B      = 499)
#' out$p.value
#' }
#' @importFrom stats glm residuals weighted.mean quantile median
#' @export
deshape_glm_resid_test <- function(formula, data,
                                   mode = c("center", "dispersion", "skewness"),
                                   family = gaussian(),
                                   group,  # grouping variable name
                                   alternative = c("two.sided", "greater", "less"),
                                   q = 0.25, B = 999, seed = NULL) {
  mode <- match.arg(mode)
  alternative <- match.arg(alternative)
  if (!is.null(seed)) set.seed(seed)
  
  # Fit GLM and extract residuals
  glm_fit <- glm(formula, data = data, family = family)
  resids <- residuals(glm_fit)
  
  g <- data[[group]]
  if (!is.factor(g)) g <- factor(g)
  G_list <- split(resids, g)
  
  # Define test statistic
  stat_fun <- switch(mode,
                     center = function(v) median(v),
                     dispersion = function(v) quantile(v, 1 - q) - quantile(v, q),
                     skewness = function(v) {
                       q_upper <- quantile(v, 1 - q)
                       q_median <- median(v)
                       q_lower <- quantile(v, q)
                       q_upper - 2 * q_median + q_lower
                     }
  )
  
  group_stats <- sapply(G_list, stat_fun)
  weights <- sapply(G_list, length)
  grand_stat <- weighted.mean(group_stats, w = weights)
  T_obs <- weighted.mean((group_stats - grand_stat)^2, w = weights)
  
  # Permutation test
  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    g_perm <- sample(g)
    G_star <- split(resids, g_perm)
    group_star_stats <- sapply(G_star, stat_fun)
    grand_star_stat <- weighted.mean(group_star_stats, w = weights)
    T_perm[b] <- weighted.mean((group_star_stats - grand_star_stat)^2, w = weights)
  }
  
  p_val <- switch(alternative,
                  "two.sided" = mean(T_perm >= T_obs),
                  "greater" = mean(T_perm >= T_obs),
                  "less" = mean(T_perm <= T_obs)
  )
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm,
    mode = mode,
    method = "GLM-residual-based DeSHAPE Test",
    B = B
  ))
}

