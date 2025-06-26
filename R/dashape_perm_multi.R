#' Permutation-Based Multi-Group Comparison
#'
#' Performs a permutation-based test across multiple groups to assess differences in
#' center (median), dispersion (IQR-based), or asymmetry (quantile-based).
#'
#' @param formula A formula of the form `response ~ group`, where the group variable has more than two levels.
#' @param data A data frame containing the variables in the formula.
#' @param mode Type of test to perform. Must be one of `"center"` (permutation median ANOVA), `"dispersion"` (IQR-based permutation test), or `"skewness"` (quantile-based asymmetry test).
#' @param perm Number of permutations. Default is 999. Used for all modes.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return This function prints the permutation p-value to the console.
#'
#' @importFrom car leveneTest
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(Shannon = rnorm(90), Group = rep(c("A", "B", "C"), each = 30))
#' deshape_perm_multi(Shannon ~ Group, data = data, mode = "center")
#' }
deshape_perm_multi <- function(formula, data, mode = c("center", "dispersion", "skewness"), perm = 999, seed) {
  mode <- match.arg(mode)
  
  vars <- all.vars(formula)
  response <- vars[1]
  predictor <- vars[2]
  
  predictor_vals <- factor(data[[predictor]])
  if (nlevels(predictor_vals) <= 2) {
    stop(paste("The variable", predictor, "should have more than two groups!"))
  }
  
  if (mode == "center") {
    result <- perm_median_anova(data, outcome = response, group = predictor, B = perm, seed = seed)
    cat("Permutation p-value:", result$p.value, "\n")
  } else if (mode == "dispersion") {
    result <- perm_dispersion_anova(data, outcome = response, group = predictor, B = perm, seed = seed)
    cat("Permutation p-value:", result$p.value, "\n")
  } else if (mode == "skewness") {
    result <- perm_asymmetry_anova(data, outcome = response, group = predictor, B = perm, seed = seed)
    cat("Permutation p-value:", result$p.value, "\n")
  } 
}

# Internal helper: not exported
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

# Internal helper: not exported
perm_dispersion_anova <- function(data, outcome, group, q = 0.25, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  x <- data[[outcome]]
  g <- data[[group]]
  
  if (!is.numeric(x)) stop("Outcome must be numeric.")
  if (!is.factor(g)) g <- factor(g)
  
  # Define dispersion (IQR) function
  dispersion <- function(v) {
    quantile(v, 1 - q, names = FALSE, type = 7) - quantile(v, q, names = FALSE, type = 7)
  }
  
  K <- nlevels(g)
  G_list <- split(x, g)
  group_disp <- sapply(G_list, dispersion)
  grand_disp <- mean(group_disp)
  T_obs <- mean((group_disp - grand_disp)^2)
  
  # Permutation step
  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    g_perm <- sample(g)
    G_star <- split(x, g_perm)
    group_star_disp <- sapply(G_star, dispersion)
    grand_star_disp <- mean(group_star_disp)
    T_perm[b] <- mean((group_star_disp - grand_star_disp)^2)
  }
  
  p_val <- mean(T_perm >= T_obs)
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm,
    method = "Permutation-based Dispersion ANOVA",
    B = B
  ))
}

# Internal helper: not exported
perm_asymmetry_anova <- function(data, outcome, group, q = 0.1, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  x <- data[[outcome]]
  g <- data[[group]]
  
  if (!is.numeric(x)) stop("Outcome must be numeric.")
  if (!is.factor(g)) g <- factor(g)
  
  # Define asymmetry score function
  asymmetry <- function(v) {
    q_upper <- quantile(v, 1 - q, names = FALSE, type = 7)
    q_median <- quantile(v, 0.5, names = FALSE, type = 7)
    q_lower <- quantile(v, q, names = FALSE, type = 7)
    return(q_upper - 2 * q_median + q_lower)
  }
  
  K <- nlevels(g)
  G_list <- split(x, g)
  group_asym <- sapply(G_list, asymmetry)
  grand_asym <- mean(group_asym)
  T_obs <- mean((group_asym - grand_asym)^2)
  
  # Permutation
  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    g_perm <- sample(g)
    G_star <- split(x, g_perm)
    group_star_asym <- sapply(G_star, asymmetry)
    grand_star_asym <- mean(group_star_asym)
    T_perm[b] <- mean((group_star_asym - grand_star_asym)^2)
  }
  
  p_val <- mean(T_perm >= T_obs)
  
  return(list(
    statistic = T_obs,
    p.value = p_val,
    T_perm = T_perm,
    method = "Permutation-based Asymmetry ANOVA",
    B = B
  ))
}
