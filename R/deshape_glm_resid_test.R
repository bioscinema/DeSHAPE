#' DeSHAPE GLM-Residual Group Shape Test
#'
#' Applies the **DeSHAPE** framework to *response* residuals from a
#' generalized-linear model (GLM) in order to compare distributional
#' *shape* features—location, dispersion, or skewness—across
#' user-defined groups.  The procedure
#' \enumerate{
#'   \item fits the GLM specified by \code{formula} and extracts raw
#'         response residuals;
#'   \item computes a group-specific summary statistic targeting the
#'         selected shape aspect (\code{"center"}, \code{"dispersion"},
#'         or \code{"skewness"});
#'   \item forms a cross-group sum-of-squares contrast statistic,
#'         \eqn{T};
#'   \item constructs a permutation reference distribution for
#'         \eqn{T} by re-shuffling group labels \code{B} times;
#'   \item returns the observed statistic together with its permutation
#'         \eqn{p}\,-value.
#' }
#' The test generalises the original DeSHAPE methodology to any GLM
#' handled by \code{\link[stats]{glm}} (e.g.\ Gaussian, Poisson,
#' binomial) and is particularly useful when covariates must be adjusted
#' prior to shape comparison.
#'
#' @section Test statistics:
#' \describe{
#'   \item{\code{mode = "center"}}{Group median of residuals.}
#'   \item{\code{mode = "dispersion"}}{Inter-quantile range
#'         \eqn{Q_{1-q} - Q_q}.}
#'   \item{\code{mode = "skewness"}}{Quantile-based third-order contrast
#'         \eqn{Q_{1-q} - 2Q_{0.5} + Q_q}.}
#' }
#'
#' @param formula A model \code{\link[stats]{formula}} passed to
#'   \code{\link[stats]{glm}}.
#' @param data    Data frame containing variables referenced in
#'   \code{formula} and the grouping variable.
#' @param mode    Character string indicating which aspect of residual
#'   shape to test; one of \code{"center"}, \code{"dispersion"},
#'   \code{"skewness"}.  See *Test statistics*.
#' @param family  A GLM family object (default
#'   \code{\link[stats]{gaussian}}).
#' @param group   **Character string** giving the column name in
#'   \code{data} that defines the groups whose residual shapes are
#'   compared.
#' @param alternative Alternative hypothesis:
#'   \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param q       Lower-tail probability used in the inter-quantile
#'   contrasts; must satisfy \eqn{0 < q < 0.5}.  Default \code{0.25}.
#' @param B       Integer; number of label permutations (default
#'   \code{999}).
#' @param seed    Optional integer seed for reproducibility.
#'
#' @return A list with components
#' \describe{
#'   \item{\code{statistic}}{Observed test statistic \eqn{T}.}
#'   \item{\code{p.value}}{Permutation \eqn{p}\,-value.}
#'   \item{\code{T_perm}}{Vector of permuted statistics
#'         (length \code{B}).}
#'   \item{\code{mode}}{Requested shape aspect.}
#'   \item{\code{method}}{Description of the test.}
#'   \item{\code{B}}{Number of permutations performed.}
#' }
#'
#' @details
#' The sum-of-squares statistic aggregates deviations of each
#' group-level summary from the overall summary.  Under the null
#' hypothesis of equal residual shape across groups, label permutation
#' yields an exact test when residuals are exchangeable and an
#' asymptotically valid test otherwise.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n  <- 120
#' g  <- gl(3, n/3, labels = c("A", "B", "C"))
#' x1 <- rnorm(n)
#' y  <- 0.5 + 1.2 * x1 + rnorm(n, sd = 0.8)      # Gaussian GLM
#' dat <- data.frame(y, x1, g)
#'
#' # Test equality of residual dispersion across groups
#' out <- deshape_glm_resid_test(y ~ x1,
#'                               data  = dat,
#'                               group  = "g",
#'                               alternative = "two.sided",
#'                               mode   = "dispersion",
#'                               B      = 499)
#' out$p.value
#' }
#' @importFrom stats glm residuals quantile median
#' @export
deshape_glm_resid_test <- function(
    formula, data,
    mode = c("center","dispersion","skewness"),
    group,
    alternative = c("two.sided","greater","less"),
    q = 0.25,
    B = 999,
    seed = NULL,
    family = gamlss.dist::JSU()   # <-- more stable default
){
  mode <- match.arg(mode)
  alternative <- match.arg(alternative)
  if(!is.null(seed)) set.seed(seed)
  
  # ---------------------------
  # 1. Remove group from formula (only regress on Z)
  # ---------------------------
  rhs <- attr(terms(formula), "term.labels")
  rhs <- rhs[rhs != group]
  if(length(rhs) == 0)
    stop("No predictors remain after dropping the group term.")
  
  lhs <- deparse(formula[[2]])
  formula_nog <- as.formula(
    paste(lhs, "~", paste(rhs, collapse = " + "))
  )
  
  # ---------------------------
  # 2. Fit GAMLSS: F(X|Z)
  # ---------------------------
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })
  
  gam_fit <- gamlss(
    formula = formula_nog,
    data = data,
    family = family,
    trace = FALSE,
    n.cyc = 200
  )
  
  if(!gam_fit$converged)
    warning("GAMLSS did not fully converge — consider changing family.")
  
  # ---------------------------
  # 3. PIT residuals U = F_hat(X | Z)
  # ---------------------------
  # pick correct p-distribution depending on family
  fam_name <- family$family[1]
  
  pfun <- switch(
    fam_name,
    "JSU"   = gamlss.dist::pJSU,
    "ST5"   = gamlss.dist::pST5,
    "SEP4"  = gamlss.dist::pSEP4,
    "NO"    = pnorm,
    stop(paste("p-function not implemented for", fam_name))
  )
  
  mu    <- fitted(gam_fit, "mu")
  sigma <- fitted(gam_fit, "sigma")
  nu    <- fitted(gam_fit, "nu")
  tau   <- fitted(gam_fit, "tau")
  
  Xv <- data[[lhs]]
  
  U <- pfun(
    q = Xv,
    mu = mu,
    sigma = sigma,
    nu = nu,
    tau = tau
  )
  
  # avoid 0/1 boundaries
  eps <- 1e-10
  U <- pmin(pmax(U, eps), 1 - eps)
  
  # ---------------------------
  # 4. Convert to Normal residual
  # ---------------------------
  R <- qnorm(U)
  
  # ---------------------------
  # 5. Compute group statistics
  # ---------------------------
  g <- data[[group]]
  if(!is.factor(g)) g <- factor(g)
  G_list <- split(R, g)
  
  stat_fun <- switch(
    mode,
    center = function(v) median(v),
    dispersion = function(v) quantile(v, 1-q) - quantile(v, q),
    skewness = function(v){
      qu <- quantile(v, 1-q)
      qm <- median(v)
      ql <- quantile(v, q)
      qu - 2*qm + ql
    }
  )
  
  group_stats <- sapply(G_list, stat_fun)
  grand_stat <- stat_fun(R)
  T_obs <- sum((group_stats - grand_stat)^2)
  
  # ---------------------------
  # 6. Permutation test (permute group labels)
  # ---------------------------
  T_perm <- numeric(B)
  for(b in seq_len(B)){
    g_perm <- sample(g)
    G_star <- split(R, g_perm)
    group_star_stats <- sapply(G_star, stat_fun)
    T_perm[b] <- sum((group_star_stats - grand_stat)^2)
  }
  
  p_val <- switch(
    alternative,
    "two.sided" = mean(abs(T_perm) >= abs(T_obs)),
    "greater"   = mean(T_perm >= T_obs),
    "less"      = mean(T_perm <= T_obs)
  )
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm,
    mode = mode,
    family = fam_name,
    method = "GAMLSS-PIT residual-based DeSHAPE Test",
    B = B
  ))
}

