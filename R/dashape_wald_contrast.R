#' Wald Contrast Test for Quantile-Regression Coefficients
#'
#' Performs a Wald-type test on linear contrasts of quantile-regression
#' coefficients across multiple quantile levels. Coefficients are estimated
#' at each \eqn{\tau \in \texttt{taus}} via \code{quantreg::rq}, stacked
#' (by \eqn{\tau}), and tested with an asymptotic covariance built from a
#' kernel density estimate of residuals at 0.
#'
#' @param formula A model formula of the form \code{response ~ predictors}.
#'   The \strong{first non-intercept term} must be a binary group variable
#'   (so that \code{model.matrix} produces exactly one column for it) when
#'   \code{mode != "customize"}.
#' @param data A data frame containing variables in \code{formula}.
#' @param mode Contrast-construction mode. One of
#'   \code{"customize"}, \code{"dispersion"}, or \code{"symmetry"}.
#'   If \code{mode = "customize"}, you must supply \code{taus} and \code{contrast}.
#'   If \code{mode = "dispersion"} or \code{"symmetry"}, \code{contrast} is
#'   constructed automatically and any user-supplied \code{contrast} is ignored.
#' @param taus Numeric vector of quantile levels in ascending order.
#'   Defaults depend on \code{mode}:
#'   \itemize{
#'     \item \code{mode = "dispersion"}: \code{c(0.25, 0.75)} (requires exactly 2 quantiles)
#'     \item \code{mode = "symmetry"}: \code{c(0.1, 0.5, 0.9)} (requires exactly 3 quantiles)
#'     \item \code{mode = "customize"}: must be supplied by the user
#'   }
#' @param contrast A numeric vector (or single-column matrix) specifying the
#'   linear contrast of the stacked coefficients \eqn{\beta(\tau_1),\ldots,\beta(\tau_K)}.
#'   Its length must be \code{p * K}, where \code{p} is \code{ncol(model.matrix(...))}
#'   and \code{K = length(taus)}. Ignored unless \code{mode = "customize"}.
#' @param alternative Alternative hypothesis. One of \code{"two.sided"},
#'   \code{"greater"}, or \code{"less"}. Defaults to \code{"two.sided"}.
#'   P-values are computed from \eqn{Z = \sqrt{\chi^2_1}} accordingly.
#' @param kernel Currently only \code{"gaussian"} is supported (placeholder for compatibility).
#'
#' @details
#' \strong{Coefficient stacking.}
#' For each \eqn{\tau_k}, let the design have \code{p} columns in the order
#' \code{(Intercept, group, covariate1, ...)}. The stacked coefficient vector is
#' \code{(β(τ1), β(τ2), ..., β(τK))} of length \code{p * K}, i.e.
#' \code{[Int τ1, group τ1, cov1 τ1, ..., Int τK, group τK, cov1 τK, ...]}.
#'
#' \strong{Auto-constructed contrasts (when \code{mode != "customize"}).}
#' Let the \emph{group} column be the first non-intercept column in \code{model.matrix}.
#' Other columns (intercept and covariates) receive weight 0 by default.
#' \itemize{
#'   \item \emph{Dispersion test} (\code{mode = "dispersion"}, \code{taus = (τ_L, τ_U)}):
#'         places \code{-1} on the group at \eqn{τ_L} and \code{+1} on the group at \eqn{τ_U}.
#'         This tests whether the group effect increases from lower to upper quantile.
#'   \item \emph{Symmetry test} (\code{mode = "symmetry"}, \code{taus = (τ_L, τ_M, τ_U)}):
#'         places \code{+1, -2, +1} on the group at \eqn{τ_L, τ_M, τ_U}, respectively,
#'         testing \eqn{\beta_{\text{group}}(τ_L) - 2\beta_{\text{group}}(τ_M) + \beta_{\text{group}}(τ_U) = 0}.
#' }
#'
#' \strong{Covariance.}
#' The covariance of the stacked coefficients is assembled using
#' \eqn{\min(\tau_k, \tau_\ell) - \tau_k \tau_\ell} blocks and
#' per-\eqn{\tau} density estimates at 0 from residuals
#' (Gaussian kernel; Silverman's rule for bandwidth).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{test_stat}: Wald chi-squared statistic with 1 degree of freedom.
#'   \item \code{p_value}: P-value computed from \eqn{Z = \sqrt{\chi^2_1}} under the
#'         specified \code{alternative}.
#' }
#'
#' @section Requirements for non-custom modes:
#' \itemize{
#'   \item The first non-intercept term in \code{formula} must be a binary group that
#'         produces exactly one column in \code{model.matrix}.
#'   \item \code{taus} must be strictly ascending; \code{length(taus)} must be 2
#'         for \code{"dispersion"} and 3 for \code{"symmetry"} (defaults provided).
#' }
#'
#' @examples
#' \dontrun{
#' library(quantreg)
#' set.seed(2025)
#' n <- 400
#' group <- rbinom(n, 1, 0.5)
#' Depth <- rnorm(n, mean = 0.7 * group, sd = 1)
#' eps0 <- rnorm(n, 0, 1)
#' eps1 <- rnorm(n, 0, 1.8) + 0.3*(rexp(n, 1) - 1)
#' y <- 0.5 * Depth + ifelse(group == 1, eps1, eps0)
#' dat <- data.frame(
#'   Shannon = y,
#'   group_prefix = factor(group, levels = c(0,1)),
#'   Depth = Depth
#' )
#'
#' ## Dispersion (auto taus/contrast): tests increase from τ=0.25 to τ=0.75
#' deshape_wald_contrast(
#'   Shannon ~ group_prefix + Depth,
#'   data = dat,
#'   mode = "dispersion"
#' )
#'
#' ## Symmetry (auto taus/contrast): tests +1·β_g(τ_L) - 2·β_g(τ_M) + 1·β_g(τ_U) = 0
#' deshape_wald_contrast(
#'   Shannon ~ group_prefix + Depth,
#'   data = dat,
#'   mode = "symmetry"
#' )
#'
#' ## Customize: provide taus and contrast explicitly (matches dispersion auto)
#' taus <- c(0.25, 0.75)
#' # p = 3 (Intercept, group_prefix, Depth); K = 2; length = 6
#' contrast <- c(0, -1, 0, 0, +1, 0)
#' deshape_wald_contrast(
#'   Shannon ~ group_prefix + Depth,
#'   data = dat,
#'   mode = "customize",
#'   taus = taus,
#'   contrast = contrast,
#'   alternative = "greater"
#' )
#' }
#'
#' @importFrom quantreg rq
#' @importFrom stats model.matrix model.response model.frame pnorm sd dnorm
#' @export
deshape_wald_contrast <- function(formula, data, 
                                  mode = c("customize","dispersion","symmetry"),q = 0.1,
                                  taus = NULL, contrast = NULL, alternative = "two.sided",
                                  kernel = "gaussian") {
  if (!all(alternative %in% c("two.sided", "greater", "less"))) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  mode <- match.arg(mode)
  if (is.null(taus)) {
    if (mode == "dispersion") {
      taus <- c(0.25, 0.75)
    } else if (mode == "symmetry") {
      taus <- c(q, 0.5, 1-q)
    }
  }
  if (mode == "customize" && is.null(taus))
    stop("When mode = 'customize', you must supply 'taus'.")
  n <- nrow(data)
  K <- length(taus)
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  p <- ncol(X)
  mf0     <- model.frame(formula, data)
  x_name  <- attr(terms(mf0), "term.labels")[1L]
  x0      <- mf0[[x_name]]
  n_group <- if (is.factor(x0)) nlevels(droplevels(x0)) else length(unique(x0))
  if (mode == "dispersion" && K != 2L)
    stop("Dispersion mode requires exactly 2 quantiles (low, high).")
  if (mode == "symmetry" && K != 3L)
    stop("Symmetry mode requires exactly 3 quantiles (low, mid, high).")
  if ((mode == "dispersion" || mode == "symmetry") && is.unsorted(taus))
    stop("'taus' must be in ascending order.")
  
  assign_vec <- attr(X, "assign")
  idx_group <- which(assign_vec == 1L)
  #if (mode != "customize" && length(idx_group) != 1L)
  #  stop("For non-custom modes, the first non-intercept term must map to exactly one design column (binary group).")
  
  if (mode != "customize" && !is.null(contrast))
    warning("'contrast' is ignored unless mode = 'customize'; auto-constructing from 'mode'.")
  
  beta_mat <- matrix(NA, nrow = p, ncol = K)
  H_list <- vector("list", K)
  S <- crossprod(X) / n  # fixed design
  
  for (k in seq_along(taus)) {
    tau <- taus[k]
    fit <- rq(y ~ X - 1, tau = tau)  # "-1" to use X as full design
    beta_hat <- coef(fit)
    beta_mat[, k] <- beta_hat
    
    resid <- y - X %*% beta_hat
    
    # Kernel density estimate at 0 using Silverman's rule
    h <- 1.06 * sd(resid) * n^(-1/5)
    f_hat <- mean(dnorm(resid / h)) / h
    
    H_list[[k]] <- f_hat * crossprod(X) / n
  }
  
  beta_vec <- as.vector(beta_mat)
  Sigma <- matrix(0, nrow = p * K, ncol = p * K)
  
  for (k in 1:K) {
    for (l in 1:K) {
      factor <- min(taus[k], taus[l]) - taus[k] * taus[l]
      Hk_inv <- solve(H_list[[k]])
      Hl_inv <- solve(H_list[[l]])
      block <- factor * Hk_inv %*% S %*% Hl_inv / n
      row_idx <- ((k - 1) * p + 1):(k * p)
      col_idx <- ((l - 1) * p + 1):(l * p)
      Sigma[row_idx, col_idx] <- block
    }
  }
  
  # ---- Build or validate contrast ------------------------------------------
  if (mode != "customize") {
    if(n_group > 2){
      contrast <- numeric(p * K)
      for (k in seq_len(K)) {
        off <- (k - 1L) * p
        if (mode == "dispersion") {
          if (length(idx_group) != 2L)
            stop("Dispersion mode expects exactly two proxy columns in the first non-intercept term.")
          if (k == 1L) contrast[off + idx_group] <- c(-1, +1)  # low tau
          if (k == 2L) contrast[off + idx_group] <- c(+1, -1)  # high tau
        } else if (mode == "symmetry") {
          if (length(idx_group) != 2L)
            stop("Symmetry mode expects exactly two proxy columns in the first non-intercept term.")
          if (k == 1L) contrast[off + idx_group] <- c(+1, -1)
          if (k == 2L) contrast[off + idx_group] <- c(-2, +2)
          if (k == 3L) contrast[off + idx_group] <- c(+1, -1)
        }
      }
    } else {
      if (length(idx_group) != 1L)
        stop("For n_group = 2, the first non-intercept term must map to exactly one design column.")
      contrast <- numeric(p * K)
      g <- idx_group[1L]
      if (mode == "dispersion") {
        contrast[g]         <- -1
        contrast[p + g]     <- +1
      } else if (mode == "symmetry") {
        contrast[g]         <- +1
        contrast[p + g]     <- -2
        contrast[2L * p + g] <- +1
      }
    }
  } else {
    if (is.null(contrast))
      stop("When mode = 'customize', you must supply 'contrast'.")
  }
  
  if (length(contrast) != p * K)
    stop("contrast must have length p * K; got ", length(contrast),
         " while p * K = ", p * K, ".")
  
  # Wald statistic
  contrast <- matrix(contrast, ncol = 1)
  est <- as.numeric(t(contrast) %*% beta_vec)
  var <- as.numeric(t(contrast) %*% Sigma %*% contrast)
  stat <- (est^2) / var
  z <- sqrt(stat)
  
  pval <- switch(alternative,
                 "two.sided" = 2 * (1 - pnorm(abs(z))),
                 "greater" = 1 - pnorm(z),
                 "less" = pnorm(z))
  
  return(list(
    test_stat = stat,
    p_value = pval
  ))
}
