#' Wald Contrast Test for Quantile Regression Coefficients
#'
#' Performs a Wald-type test on contrasts of estimated quantile regression coefficients
#' across multiple quantile levels. The asymptotic covariance is estimated using kernel
#' density estimation of residuals.
#'
#' @param formula A formula of the form `response ~ predictors`.
#' @param data A data frame containing the variables in the model.
#' @param taus A numeric vector of quantile levels (e.g., `c(0.25, 0.5, 0.75)`).
#' @param contrast A contrast vector or matrix specifying linear combinations of the stacked coefficients across quantiles.
#' @param alternative A string specifying the alternative hypothesis. One of `"two.sided"`, `"greater"`, or `"less"`. Default is `"two.sided"`.
#' @param kernel Currently only `"gaussian"` is supported. Included for compatibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{test_stat}: The Wald chi-squared test statistic.
#'   \item \code{p_value}: The p-value based on the normal approximation of the Wald statistic.
#' }
#'
#' @importFrom quantreg rq
#' @export
#'
#' @examples
#' \dontrun{
#' library(quantreg)
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n)
#' y <- 1 + 2 * x + rnorm(n)
#' df <- data.frame(y = y, x = x)
#' taus <- c(0.25, 0.5, 0.75)
#' contrast <- matrix(c(0, 1, 0, -1, 0, 0), ncol = 1)  # Compare slopes at tau=0.25 and 0.5
#' deshape_wald_contrast(y ~ x, data = df, taus = taus, contrast = contrast)
#' }
deshape_wald_contrast <- function(formula, data, 
                                  mode = c("customize","dispersion","symmetry"),
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
      taus <- c(0.1, 0.5, 0.9)
    }
  }
  if (mode == "customize" && is.null(taus))
    stop("When mode = 'customize', you must supply 'taus'.")
  n <- nrow(data)
  K <- length(taus)
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  p <- ncol(X)
  
  if (mode == "dispersion" && K != 2L)
    stop("Dispersion mode requires exactly 2 quantiles (low, high).")
  if (mode == "symmetry" && K != 3L)
    stop("Symmetry mode requires exactly 3 quantiles (low, mid, high).")
  if ((mode == "dispersion" || mode == "symmetry") && is.unsorted(taus))
    stop("'taus' must be in ascending order.")
  
  assign_vec <- attr(X, "assign")
  idx_group <- which(assign_vec == 1L)
  if (mode != "customize" && length(idx_group) != 1L)
    stop("For non-custom modes, the first non-intercept term must map to exactly one design column (binary group).")
  
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
    # Auto-construct contrast that targets the first non-intercept (binary group) column
    contrast <- numeric(p * K)
    for (k in seq_len(K)) {
      off <- (k - 1L) * p
      if (mode == "dispersion") {
        # K == 2: lower tau gets -1 on group; upper tau gets +1
        if (k == 1L) contrast[off + idx_group] <- -1
        if (k == K)  contrast[off + idx_group] <- +1
      } else if (mode == "symmetry") {
        # K == 3: weights +1, -2, +1 across taus on the group column
        if (k == 1L) contrast[off + idx_group] <- +1
        if (k == 2L) contrast[off + idx_group] <- -2
        if (k == 3L) contrast[off + idx_group] <- +1
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
