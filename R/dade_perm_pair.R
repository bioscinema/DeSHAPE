#' Permutation-Based Group Comparison (Binary)
#'
#' Performs a permutation-based statistical test to compare two groups on center (median),
#' dispersion (IQR-based Brown-Forsythe variant), or asymmetry (quantile skewness).
#'
#' @param formula A formula of the form `response ~ group`, where the group variable has exactly two levels.
#' @param data A data frame containing the variables in the formula.
#' @param mode Type of test to perform. Must be one of `"center"`, `"dispersion"`, or `"skewness"`.
#' @param alternative A character string specifying the alternative hypothesis.
#'   One of `"greater"`, `"less"`, or `"two.sided"`. Only used for the center test.
#' @param perm Number of permutations. Default is 999.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return This function prints the permutation p-value to the console.
#'
#' @importFrom car leveneTest
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(Shannon = rnorm(100), Group = rep(c("A", "B"), 50))
#' dade_perm_pair(Shannon ~ Group, data = data, mode = "center")
#' }
dade_perm_pair <- function(formula, data, 
                           mode = c("center", "dispersion", "skewness"), 
                           alternative = c("greater", "less", "two.sided"), 
                           perm = 999, seed = NULL) {
  
  mode <- match.arg(mode)
  alternative <- match.arg(alternative)
  
  vars <- all.vars(formula)
  response <- vars[1]
  predictor <- vars[2]
  
  predictor_vals <- factor(data[[predictor]])
  if (nlevels(predictor_vals) != 2) {
    stop(paste("The variable", predictor, "should have exactly two groups!"))
  }
  
  levels_group <- levels(predictor_vals)
  x <- data[[response]][data[[predictor]] == levels_group[1]]
  y <- data[[response]][data[[predictor]] == levels_group[2]]
  
  if (mode == "center") {
    result <- perm_median_test(x = x, y = y, B = perm, alternative = alternative, seed = seed)
    cat("Permutation p-value:", result$p.value, "\n")
    
  } else if (mode == "dispersion") {
    result <- perm_dispersion_test(x = x, y = y, B = perm, alternative = alternative, seed = seed)
    cat("Permutation p-value:", result$p.value, "\n")
  } else if (mode == "skewness") {
    result <- perm_asymmetry_test(x = x, y = y, B = perm, alternative = alternative, seed = seed)
    cat("Permutation p-value:", result$p.value, "\n")
  } 
}

# Internal helper: not exported
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

# Internal helper: not exported
perm_dispersion_test <- function(x, y, q = 0.25, B = 1000, seed = NULL,
                                 alternative = c("two.sided", "greater", "less")) {
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors")
  if (q <= 0 || q >= 0.5) stop("q must be in (0, 0.5)")
  alternative <- match.arg(alternative)
  if (!is.null(seed)) set.seed(seed)
  
  # Step 1: Define dispersion function (IQR at level q)
  iqr_score <- function(v, q = 0.25) {
    q_upper <- quantile(v, 1 - q, names = FALSE, type = 7)
    q_lower <- quantile(v, q, names = FALSE, type = 7)
    return(q_upper - q_lower)
  }
  
  D_x <- iqr_score(x)
  D_y <- iqr_score(y)
  T_obs <- D_x - D_y  
  
  # Step 2: Permutation
  Z <- c(x, y)
  n_x <- length(x)
  T_perm <- numeric(B)
  
  for (b in seq_len(B)) {
    permuted <- sample(Z)
    x_star <- permuted[1:n_x]
    y_star <- permuted[(n_x + 1):length(Z)]
    D_x_star <- iqr_score(x_star)
    D_y_star <- iqr_score(y_star)
    T_perm[b] <- D_x_star - D_y_star
  }
  
  # Step 3: Compute p-value
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
    T_perm = T_perm,
    alternative = alternative
  ))
}

# Internal helper: not exported
perm_asymmetry_test <- function(x, y, q = 0.05, B = 1000, seed = NULL,
                                alternative = c("two.sided", "greater", "less")) {
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors")
  if (q <= 0 || q >= 0.5) stop("q must be in (0, 0.5)")
  alternative <- match.arg(alternative)
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
  T_obs <- A_x - A_y  
  
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
    T_perm[b] <- A_x_star - A_y_star
  }
  
  # Step 3: Compute p-value based on alternative
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
    T_perm = T_perm,
    alternative = alternative
  ))
}
