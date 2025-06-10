#' Permutation-Based Multi-Group Comparison
#'
#' Performs a permutation-based test across multiple groups to assess differences in
#' center (median) or dispersion (median-based Levene test).
#'
#' @param formula A formula of the form `response ~ group`, where the group variable has more than two levels.
#' @param data A data frame containing the variables in the formula.
#' @param mode Type of test to perform. Must be one of `"center"` (permutation median ANOVA) or `"dispersion"` (Levene test).
#' @param perm Number of permutations. Default is 999. Used only for `mode = "center"`.
#'
#' @return The function returns:
#' \itemize{
#'   \item If `mode = "center"`: A list with the observed test statistic, p-value, and permutation distribution.
#'   \item If `mode = "dispersion"`: A `data.frame` with the Levene test results from `car::leveneTest`.
#' }
#'
#' @importFrom car leveneTest
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(Shannon = rnorm(90), Group = rep(c("A", "B", "C"), each = 30))
#' dade_perm_multi(Shannon ~ Group, data = data, mode = "center")
#' }
dade_perm_multi <- function(formula, data, mode = c("center", "dispersion"), perm = 999) {
  mode <- match.arg(mode)
  
  vars <- all.vars(formula)
  response <- vars[1]
  predictor <- vars[2]
  
  predictor_vals <- factor(data[[predictor]])
  if (nlevels(predictor_vals) <= 2) {
    stop(paste("The variable", predictor, "should have more than two groups!"))
  }
  
  if (mode == "center") {
    result <- perm_median_anova(data, outcome = response, group = predictor, B = perm, seed = 2025)
    cat("Permutation p-value:", result$p.value, "\n")
    
  } else if (mode == "dispersion") {
    data[[predictor]] <- factor(data[[predictor]])
    result <- car::leveneTest(formula, data = data, center = median)
    return(result)
  }
}

# Internal helper: not exported
perm_median_anova <- function(data, outcome, group, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  x <- data[[outcome]]
  g <- data[[group]]
  
  if (!is.numeric(x)) stop("Outcome must be numeric.")
  if (!is.factor(g)) g <- factor(g)
  
  G_list <- split(x, g)
  group_medians <- sapply(G_list, median)
  grand_median <- mean(group_medians)
  
  T_obs <- mean((group_medians - grand_median)^2)
  
  T_perm <- numeric(B)
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
