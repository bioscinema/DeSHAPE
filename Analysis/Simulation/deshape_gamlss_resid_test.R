deshape_gamlss_resid_test <- function(
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

