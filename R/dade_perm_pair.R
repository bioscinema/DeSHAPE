#' Permutation-Based Group Comparison (Binary)
#'
#' Performs a permutation-based statistical test to compare two groups on center (median),
#' dispersion (Brown-Forsythe test), or asymmetry (quantile skewness).
#'
#' @param formula A formula of the form `response ~ group`, where the group variable has exactly two levels.
#' @param data A data frame containing the variables in the formula.
#' @param mode Type of test to perform. Must be one of `"center"`, `"dispersion"`, or `"skewness"`.
#' @param alternative A character string specifying the alternative hypothesis.
#'   One of `"greater"`, `"less"`, or `"two.sided"`. Only used for the center test.
#' @param perm Number of permutations. Default is 999.
#'
#' @return The function returns test results depending on the selected mode:
#' \itemize{
#'   \item If `mode = "center"`: A list with statistic and p-value from the permutation median test.
#'   \item If `mode = "dispersion"`: A `data.frame` from `car::leveneTest`.
#'   \item If `mode = "skewness"`: A list with statistic and p-value from the asymmetry test.
#' }
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
                           perm = 999) {
  
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
    result <- perm_median_test(x, y, B = perm, alternative = alternative)
    cat("Observed median difference(", levels_group[1],"-",levels_group[2],"):", result$statistic, "\n")
    cat("Permutation p-value:", result$p.value, "\n")
    
  } else if (mode == "dispersion") {
    data[[predictor]] <- factor(data[[predictor]])
    result <- car::leveneTest(formula, data = data, center = median)
    return(result)
    
  } else if (mode == "skewness") {
    result <- perm_asymmetry_test(x, y, q = 0.05, B = perm)
    cat("Asymmetry difference:", result$statistic, "\n")
    cat("Permutation p-value:", result$p.value, "\n")
  }
}
