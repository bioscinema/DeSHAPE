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
#'   \item \code{test_stat}: The Wald test statistic.
#'   \item \code{p_value}: The corresponding p-value based on normal approximation.
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
#' wald_contrast_test(y ~ x, df, taus = taus, contrast = contrast)
#' }
wald_contrast_test <- function(formula, data, taus,
                               contrast, alternative = "two.sided",
                               kernel = "gaussian") {
  if (!all(alternative %in% c("two.sided", "greater", "less"))) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  n <- nrow(data)
  K <- length(taus)
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  p <- ncol(X)
  
  beta_mat <- matrix(NA, nrow = p, ncol = K)
  H_list <- vector("list", K)
  S <- crossprod(X) / n  # fixed design
  
  for (k in seq_along(taus)) {
    tau <- taus[k]
    fit <- quantreg::rq(y ~ X - 1, tau = tau)  # use X directly
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
  
  contrast <- as.matrix(contrast)
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
