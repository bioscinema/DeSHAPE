.libPaths(c("/home/zhang.16383/RPackage", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
rep_id <- as.numeric(args[1]); if (is.na(rep_id)) rep_id <- 1L
if (rep_id < 1L) stop("rep_id must be >= 1")

# -------- Helper function
base_dir <- "result"
subdirs  <- c("sim_summary_median", "sim_summary_median_confounder", "sim_summary_sd", "sim_summary_sd_confounder",
              "sim_summary_skew", "sim_summary_skew_confounder")

if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

for (d in file.path(base_dir, subdirs)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


library(moments)
library(DeSHAPE)

SIRGamma <- function(n, a, alpha, gamma,
                     q_probs = c(0.1, 0.7, 0.9),
                     target_prop = c(0.1, 0.2, 0.1),
                     n_pool = 1e5, seed = NULL) {
  #-------------------------------------------------------------
  # Importance resampling from location-shifted Gamma
  # X = a + Gamma(alpha, gamma)
  # Ensures fixed proportions of samples in quantile intervals:
  #   <= Q(0.1) : target_prop[1]
  #   (Q(0.7), Q(0.9)] : target_prop[2]
  #   >= Q(0.9) : target_prop[3]
  #-------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  pool <- a + rgamma(n_pool, shape = alpha, scale = gamma)
  
  q_vals <- a + gamma * qgamma(q_probs, shape = alpha, scale = 1)
  names(q_vals) <- paste0("Q", q_probs)
  
  idx_low  <- which(pool <= q_vals["Q0.1"])
  idx_mid  <- which(pool > q_vals["Q0.7"] & pool <= q_vals["Q0.9"])
  idx_high <- which(pool >= q_vals["Q0.9"])
  
  n_low  <- round(target_prop[1] * n)
  n_mid  <- round(target_prop[2] * n)
  n_high <- round(target_prop[3] * n)
  n_rest <- n - (n_low + n_mid + n_high)
  
  # importance resampling
  samp_low  <- sample(pool[idx_low],  n_low,  replace = TRUE)
  samp_mid  <- sample(pool[idx_mid],  n_mid,  replace = TRUE)
  samp_high <- sample(pool[idx_high], n_high, replace = TRUE)
  
  rest_idx <- setdiff(seq_along(pool), c(idx_low, idx_mid, idx_high))
  samp_rest <- sample(pool[rest_idx], n_rest, replace = TRUE)
  
  final_sample <- sample(c(samp_low, samp_mid, samp_high, samp_rest))
  
  list(
    sample = final_sample,
    quantiles = q_vals,
    summary = data.frame(
      region = c("<=Q0.1", "(Q0.7,Q0.9]", ">=Q0.9", "rest"),
      target_prop = c(target_prop, 1 - sum(target_prop)),
      realized_prop = c(
        mean(final_sample <= q_vals["Q0.1"]),
        mean(final_sample > q_vals["Q0.7"] & final_sample <= q_vals["Q0.9"]),
        mean(final_sample >= q_vals["Q0.9"]),
        mean(final_sample > q_vals["Q0.1"] & final_sample <= q_vals["Q0.7"])
      )
    )
  )
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
### case 1: median difference (two group)

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 20, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    a = 1
    alpha = 4
    gamma = 1/2
    
    G1 = SIRGamma(n,a,alpha,gamma)
    G2 = SIRGamma(n,a*delta0,alpha,gamma)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c(G1$sample,G2$sample)
    )
    Z <- c(rnorm(n, mean = 0, sd = 1),
           rnorm(n, mean = 1, sd = 1))
    X_median_conf <- add_confounder(X_latent = simdata$X, Z, group = simdata$group,
                                    beta_m = 1, beta_s = 0, beta_k = 0)
    simdata$X <- X_median_conf
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
                                            mode = "center", alternative = "two.sided"),
                          error = function(e) NA_real_)
    p_deshape_glm <- tryCatch(deshape_glm_resid_test(X ~ group+Z, group = "group",
                                            data = simdata, mode = "center",alternative = "two.sided")$p.value,
                              error = function(e) NA_real_)
    p_ks_glm <- tryCatch(ks_glm_resid_test(X ~ group+Z, group = "group",
                                  data = simdata, mode = "center",alternative = "two.sided"),
                         error = function(e) NA_real_)
    qr_fit <- summary(quantreg::rq(X ~ group + Z,
                           data = simdata,
                           tau = 0.5,
                           method = "fn"),se = "boot")
    p_qr <- qr_fit$coefficients["group", "Pr(>|t|)"]
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
      p_qr = p_qr,
      stringsAsFactors = FALSE
    )
  }
}

sim_summary_median <- do.call(rbind, results)

sim_summary_median$rep_id <- rep_id
write.csv(
  sim_summary_median,
  file.path(base_dir, "sim_summary_median_confounder", sprintf("sim_summary_rep%03d.csv", rep_id)),
  row.names = FALSE
)


###################################################
### case 1: median difference (two group) - non confounder

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 20, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    a = 1
    alpha = 4
    gamma = 1/2
    
    G1 = SIRGamma(n,a,alpha,gamma)
    G2 = SIRGamma(n,a*delta0,alpha,gamma)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c(G1$sample,G2$sample)
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
                                            mode = "center", alternative = "two.sided"),
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

sim_summary_median <- do.call(rbind, results)

sim_summary_median$rep_id <- rep_id
write.csv(
  sim_summary_median,
  file.path(base_dir, "sim_summary_median", sprintf("sim_summary_rep%03d.csv", rep_id)),
  row.names = FALSE
)

###################################################
### case 2: dispersion difference (two group)

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 20, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    a = 1
    alpha = 4
    gamma1 = 1
    gamma2 = 1*delta0
    
    G1 = SIRGamma(n,a,alpha,gamma1)
    G2 = SIRGamma(n,a,alpha,gamma2)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c(G1$sample-mean(G1$sample),G2$sample-mean(G2$sample))
    )
    Z <- (c(rnorm(n, mean = 0, sd = 1),
                 rnorm(n, mean = 1, sd = 1)))
    X_sd_conf <- add_confounder(X_latent = simdata$X, Z, group = simdata$group,
                                beta_m = 0.0, beta_s = 1, beta_k = 0)
    simdata$X <- X_sd_conf
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
                                            mode = "dispersion", alternative = "two.sided"),
                          error = function(e) NA_real_)
    p_deshape_glm <- tryCatch(deshape_glm_resid_test(X ~ group+Z, group = "group",
                                                     data = simdata, mode = "center",alternative = "two.sided")$p.value,
                              error = function(e) NA_real_)
    p_ks_glm <- tryCatch(ks_glm_resid_test(X ~ group+Z, group = "group",
                                  data = simdata, mode = "dispersion",alternative = "two.sided"), error = function(e) NA_real_)
    p_deshape_deconfounder <- deshape_wald_contrast(X~group+Z,data = simdata,mode = "dispersion",alternative = "two.sided")$p_value
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
  file.path(base_dir, "sim_summary_sd_confounder", sprintf("sim_summary_rep%03d.csv", rep_id)),
  row.names = FALSE
)

###################################################
### case 2: dispersion difference (two group) - no confounder

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 20, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    a = 1
    alpha = 4
    gamma1 = 1
    gamma2 = 1*delta0
    
    G1 = SIRGamma(n,a,alpha,gamma1)
    G2 = SIRGamma(n,a,alpha,gamma2)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c(G1$sample-mean(G1$sample),G2$sample-mean(G2$sample))
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
                                            mode = "dispersion", alternative = "two.sided"),
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
  file.path(base_dir, "sim_summary_sd", sprintf("sim_summary_rep%03d.csv", rep_id)),
  row.names = FALSE
)


###################################################
### case 3: skewness difference (two group)

n_vals  <- seq(10, 500, by = 20)
deltas  <- seq(1, 20, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    a = 1
    alpha1 = 1/4
    alpha2 = alpha1*delta0
    gamma = 2
    
    G1 = SIRGamma(n,a,alpha1,gamma)
    G2 = SIRGamma(n,a,alpha2,gamma)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c((G1$sample-mean(G1$sample))/sd(G1$sample),(G2$sample-mean(G2$sample))/sd(G2$sample))
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
deltas  <- seq(1, 20, by = 0.1)

results <- vector("list", length(n_vals) * length(deltas))
idx <- 0L

for (n in n_vals) {
  for (delta0 in deltas) {
    idx <- idx + 1L
    
    a = 1
    alpha1 = 1/4
    alpha2 = alpha1*delta0
    gamma = 2
    
    G1 = SIRGamma(n,a,alpha1,gamma)
    G2 = SIRGamma(n,a,alpha2,gamma)
    simdata = data.frame(
      group = c(rep(1,n),rep(2,n)),
      X = c((G1$sample-mean(G1$sample))/sd(G1$sample),(G2$sample-mean(G2$sample))/sd(G2$sample))
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