.libPaths(c("/home/zhang.16383/RPackage", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
rep_id <- as.numeric(args[1]); if (is.na(rep_id)) rep_id <- 1L
if (rep_id < 1L) stop("rep_id must be >= 1")

# -------- Helper function
base_dir <- "resultskew"
subdirs  <- c("sim_summary_skew", "sim_summary_skew_confounder")

if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

for (d in file.path(base_dir, subdirs)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


library(moments)
library(DeSHAPE)

rbeta_n <- function(n, alpha, beta, varscale = 10) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 0,
            is.numeric(alpha), length(alpha) == 1, alpha > 0,
            is.numeric(beta),  length(beta)  == 1, beta  > 0)
  betadata <- stats::rbeta(n = as.integer(n), shape1 = alpha, shape2 = beta)
  standardizeddata <- (betadata - median(betadata))/sd(betadata)
  standardizeddata * varscale
}


add_confounder <- function(X_latent, Z, group,
                           beta_m = 0, beta_s = 0, beta_k = 0) {
  
  X <- as.numeric(X_latent)
  Z <- as.numeric(Z)
  g <- as.character(group)
  
  stopifnot(length(X) == length(Z), length(Z) == length(g))
  
  centers <- tapply(X, g, median, na.rm = TRUE)
  
  cvec <- as.numeric(centers[ match(g, names(centers)) ])
  
  if (beta_m != 0) {
    X <- X + beta_m * Z
    centers <- tapply(X, g, median, na.rm = TRUE)
    cvec <- as.numeric(centers[ match(g, names(centers)) ])
  }
  
  if (beta_s != 0) {
    sZ <- as.numeric(exp(beta_s * Z))  
    X  <- cvec + sZ * (X - cvec)
  }
  
  if (beta_k != 0) {
    aZ <- as.numeric(exp(beta_k * Z))
    d  <- X - cvec
    X  <- cvec + aZ * pmax(d, 0) + (1 / aZ) * pmin(d, 0)
  }
  
  X
}



ks_glm_resid_test <- function(
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
  x <- G_list[[1]]
  y <- G_list[[2]]
  
  ks_res <- ks.test(x, y, alternative = alternative)
  
  ks_res$p.value
}
# ------- Non Confounder

set.seed(1000 + 10 * rep_id)


###################################################
### case 3: skewness difference (two group)

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 10, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    
    
    G1 = rbeta_n(n = as.integer(n), alpha = 20, beta = 20 + delta0)
    G2 = rbeta_n(n = as.integer(n), alpha = 20, beta = 20 - delta0)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c(G1,G2)
    )
    Z <- (c(rnorm(n, mean = 0, sd = 1),
                 rnorm(n, mean = 1, sd = 1)))
    X_skew_conf <- add_confounder(X_latent = simdata$X, Z, group = simdata$group,
                                  beta_m = 0.0, beta_s = 0.0, beta_k = 2)
    simdata$X <- X_skew_conf
    simdata$Z <- Z
    
    # tests
    p_wilcox <- tryCatch(wilcox.test(X ~ group, data = simdata, alternative = "two.sided")$p.value,
                         error = function(e) NA_real_)
    p_ttest  <- tryCatch(t.test(X ~ group, data = simdata, alternative = "two.sided")$p.value,
                         error = function(e) NA_real_)
    p_ks     <- tryCatch(ks.test(X ~ group, data = simdata,
                                 alternative = "two.sided")$p.value,
                         error = function(e) NA_real_)
    p_deshape <- tryCatch(deshape_perm_pair(X ~ group, data = simdata,
                                            mode = "skewness", alternative = "two.sided"),
                          error = function(e) NA_real_)
    p_deshape_glm <- tryCatch(deshape_glm_resid_test(X ~ group+Z, group = "group",
                                            data = simdata, mode = "skewness",alternative = "two.sided",q=0.1)$p.value,
                              error = function(e) NA_real_)
    p_ks_glm <- tryCatch(ks_glm_resid_test(X ~ group+Z, group = "group",
                                  data = simdata, mode = "skewness",alternative = "two.sided"), error = function(e) NA_real_)
    p_deshape_deconfounder <- deshape_wald_contrast(X~group+Z,data = simdata,mode = "symmetry",alternative = "two.sided")$p_value
    
    # store row
    results[[idx]] <- data.frame(
      n1 = n, n2 = n, delta0 = delta0,
      
      # p-values
      p_wilcox = p_wilcox,
      p_ttest  = p_ttest,
      p_ks     = p_ks,
      p_deshape = p_deshape,
      p_deshape_glm = p_deshape_glm,
      p_ks_glm = p_ks_glm,
      p_deshape_deconfounder = p_deshape_deconfounder,
      stringsAsFactors = FALSE
    )
  }
}

sim_summary_sd <- do.call(rbind, results)

sim_summary_sd$rep_id <- rep_id
write.csv(
  sim_summary_sd,
  file.path(base_dir, "sim_summary_skew_confounder", sprintf("sim_summary_rep%03d.csv", rep_id)),
  row.names = FALSE
)



###################################################
### case 3: skewness difference (two group) - non confounder

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 10, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    
    
    G1 = rbeta_n(n = as.integer(n), alpha = 20, beta = 20 + delta0)
    G2 = rbeta_n(n = as.integer(n), alpha = 20, beta = 20 - delta0)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c(G1,G2)
    )
    
    # tests
    p_wilcox <- tryCatch(wilcox.test(X ~ group, data = simdata, alternative = "two.sided")$p.value,
                         error = function(e) NA_real_)
    p_ttest  <- tryCatch(t.test(X ~ group, data = simdata, alternative = "two.sided")$p.value,
                         error = function(e) NA_real_)
    p_ks     <- tryCatch(ks.test(X ~ group, data = simdata,
                                 alternative = "two.sided")$p.value,
                         error = function(e) NA_real_)
    p_deshape <- tryCatch(deshape_perm_pair(X ~ group, data = simdata,
                                            mode = "skewness", alternative = "two.sided"),
                          error = function(e) NA_real_)
    # store row
    results[[idx]] <- data.frame(
      n1 = n, n2 = n, delta0 = delta0,
      
      # p-values
      p_wilcox = p_wilcox,
      p_ttest  = p_ttest,
      p_ks     = p_ks,
      p_deshape = p_deshape,
      
      stringsAsFactors = FALSE
    )
  }
}

sim_summary_sd <- do.call(rbind, results)

sim_summary_sd$rep_id <- rep_id
write.csv(
  sim_summary_sd,
  file.path(base_dir, "sim_summary_skew", sprintf("sim_summary_rep%03d.csv", rep_id)),
  row.names = FALSE
)