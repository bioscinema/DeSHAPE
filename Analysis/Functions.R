perm_median_test <- function(x, y, B = 1000, alternative = "two.sided", seed = NULL) {
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("alternative must be one of 'two.sided', 'greater', 'less'")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  T_obs <- median(x) - median(y)
  T_perm <- numeric(B)
  Z <- c(x, y)
  n_x <- length(x)
  
  for (b in seq_len(B)) {
    perm_labels <- sample(Z)
    x_star <- perm_labels[1:n_x]
    y_star <- perm_labels[(n_x + 1):length(Z)]
    T_perm[b] <- median(x_star) - median(y_star)
  }
  
  # Compute p-value
  if (alternative == "two.sided") {
    p_val <- mean(abs(T_perm) >= abs(T_obs))
  } else if (alternative == "greater") {
    p_val <- mean(T_perm >= T_obs)
  } else {  # alternative == "less"
    p_val <- mean(T_perm <= T_obs)
  }
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm
  ))
}


perm_median_anova <- function(data, outcome, group, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  x <- data[[outcome]]
  g <- data[[group]]
  
  if (!is.numeric(x)) stop("Outcome must be numeric.")
  if (!is.factor(g)) g <- factor(g)
  
  K <- nlevels(g)
  G_list <- split(x, g)
  group_medians <- sapply(G_list, median)
  grand_median <- mean(group_medians)
  
  T_obs <- mean((group_medians - grand_median)^2)
  
  T_perm <- numeric(B)
  n <- length(x)
  for (b in seq_len(B)) {
    g_perm <- sample(g)  # permute labels
    G_star <- split(x, g_perm)
    group_star_medians <- sapply(G_star, median)
    grand_star_median <- mean(group_star_medians)
    T_perm[b] <- mean((group_star_medians - grand_star_median)^2)
  }
  
  p_val <- mean(T_perm >= T_obs)
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm,
    method = "Permutation-based median ANOVA",
    B = B
  ))
}


perm_asymmetry_test <- function(x, y, q = 0.05, B = 1000, seed = NULL) {
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors")
  if (q <= 0 || q >= 0.5) stop("q must be in (0, 0.5)")
  if (!is.null(seed)) set.seed(seed)
  
  # Step 1: Define asymmetry function
  asym_score <- function(v) {
    q_upper <- quantile(v, 1 - q, names = FALSE, type = 7)
    q_median <- quantile(v, 0.5, names = FALSE, type = 7)
    q_lower <- quantile(v, q, names = FALSE, type = 7)
    return(q_upper - 2 * q_median + q_lower)
  }
  
  A_x <- asym_score(x)
  A_y <- asym_score(y)
  T_obs <- abs(A_x - A_y)
  
  # Step 2: Permutation test
  Z <- c(x, y)
  n_x <- length(x)
  T_perm <- numeric(B)
  
  for (b in seq_len(B)) {
    permuted <- sample(Z)
    x_star <- permuted[1:n_x]
    y_star <- permuted[(n_x + 1):length(Z)]
    A_x_star <- asym_score(x_star)
    A_y_star <- asym_score(y_star)
    T_perm[b] <- abs(A_x_star - A_y_star)
  }
  
  p_val <- mean(T_perm >= T_obs)
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm
  ))
}


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
  
  # Wald statistic
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


bootstrap_contrast_test <- function(formula, data, taus,
                                    contrast, B = 1000,
                                    alternative = "two.sided",
                                    method = "resample") {
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  n <- nrow(X)
  p <- ncol(X)
  K <- length(taus)
  
  # Point estimate
  beta_mat <- sapply(taus, function(tau) coef(rq(y ~ X - 1, tau = tau)))
  beta_vec <- as.vector(beta_mat)
  contrast <- as.matrix(contrast)
  
  delta_obs <- as.numeric(t(contrast) %*% beta_vec)
  
  # Bootstrap distribution
  delta_boot <- numeric(B)
  
  for (b in 1:B) {
    idx <- sample(1:n, replace = TRUE)
    yb <- y[idx]
    Xb <- X[idx, , drop = FALSE]
    
    beta_b <- sapply(taus, function(tau) coef(rq(yb ~ Xb - 1, tau = tau)))
    beta_b_vec <- as.vector(beta_b)
    delta_boot[b] <- as.numeric(t(contrast) %*% beta_b_vec)
  }
  
  delta_mean <- mean(delta_boot)
  delta_sd <- sd(delta_boot)
  
  # Z-like statistic
  z <- (delta_obs - delta_mean) / delta_sd
  
  pval <- switch(alternative,
                 "two.sided" = 2 * (1 - pnorm(abs(z))),
                 "greater" = 1 - pnorm(z),
                 "less" = pnorm(z)
  )
  
  return(list(
    delta_obs = delta_obs,
    bootstrap_mean = delta_mean,
    bootstrap_sd = delta_sd,
    test_stat = z^2,
    p_value = pval
  ))
}

