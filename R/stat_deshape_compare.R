#' Add DeSHAPE three-part comparison p-values to a ggplot
#'
#' @description
#' `stat_deshape_compare()` computes and displays three p-values summarizing
#' group differences in **center** (location), **dispersion**, and **asymmetry**
#' for the response mapped to `y` across the groups mapped to `x`.
#'
#' - **Unadjusted mode** (no confounders): uses your permutation engine
#'   \code{\link{deshape_perm_multi}} with \code{mode = "center"},
#'   \code{"dispersion"}, and \code{"skewness"}.
#' - **Adjusted mode** (with confounders): tests center via a median
#'   (\eqn{\tau=0.5}) quantile-regression ANOVA, and tests dispersion and
#'   asymmetry via Wald-contrast tests from
#'   \code{\link{deshape_wald_contrast}} with \code{mode = "dispersion"}
#'   and \code{"symmetry"} respectively.
#'
#' The stat draws a **single text label per panel** combining the three
#' p-values as lines (or with a custom separator), positioned using either
#' normalized parent coordinates or absolute data coordinates.
#'
#' @section Aesthetics:
#' `stat_deshape_compare()` requires the following aesthetics:
#' \itemize{
#'   \item \code{x} — grouping variable (treated categorically)
#'   \item \code{y} — response variable
#' }
#' It understands the following aesthetics (inherited by the text geom):
#' \code{hjust}, \code{vjust}, and any text aesthetics accepted by the chosen geom.
#'
#' @inheritParams ggplot2::layer
#'
#' @param confounder `NULL` or a character vector of column names giving
#'   confounders to adjust for. If supplied, you must also provide
#'   `confounder_data` containing these columns; rows will be aligned to the
#'   layer data by integer rownames if available, otherwise sequentially (see
#'   Details).
#' @param perm Integer; number of permutations used by the unadjusted
#'   permutation tests (center/dispersion/skewness).
#' @param seed Integer random seed for reproducibility.
#' @param label.sep String used to join the three lines (default `"\n"`).
#' @param label.prefix Character vector of length 3 giving the prefixes for the
#'   three reported components (defaults to `c("Center","Dispersion","Asymmetry")`).
#' @param label.x.npc, label.y.npc
#'   Position of the label expressed in normalized parent coordinates or keywords.
#'   For `x`: one of `"left"`, `"right"`, `"center"` or a numeric in `[0,1]`.
#'   For `y`: one of `"top"`, `"bottom"`, `"center"` or a numeric in `[0,1]`.
#'   Ignored if `label.x` / `label.y` are provided.
#' @param label.x, label.y
#'   Absolute data coordinates at which to draw the label. If provided, they
#'   override `label.x.npc` / `label.y.npc`.
#' @param vjust Numeric offset added to the computed vertical justification for
#'   fine-tuning label placement (useful with \code{geom = "text"}).
#' @param confounder_data A data frame providing the columns named in
#'   `confounder`. Must have at least as many rows as the panel data. See Details
#'   for row alignment rules.
#' @param ... Additional parameters passed to the chosen `geom` (e.g., size,
#'   fontface, alpha).
#'
#' @details
#' \strong{Modeling and tests.}
#' \itemize{
#'   \item Unadjusted tests call \code{\link{deshape_perm_multi}} three times
#'   on the model \code{y ~ x} with \code{mode = "center"}, \code{"dispersion"},
#'   and \code{"skewness"}.
#'   \item Adjusted tests build \code{y ~ x + confounders}. Center is tested by
#'   comparing a median-quantile regression (\eqn{\tau=0.5}) full model against
#'   a null model with confounders only via \code{stats::anova} on
#'   \code{quantreg::rq} fits with \code{joint = TRUE}. Dispersion and asymmetry
#'   are tested using \code{\link{deshape_wald_contrast}} with
#'   \code{mode = "dispersion"} and \code{mode = "symmetry"}.
#' }
#'
#' \strong{Row alignment for \code{confounder_data}.}
#' If the panel data retain integer rownames corresponding to the original input
#' rows, those indices are used to align with `confounder_data`. Otherwise, a
#' sequential fallback is used (first \eqn{n} rows), which may be inaccurate
#' under faceting/filtering. For robust alignment across facets, include a
#' stable row id in your dataset and ensure it is preserved into the layer.
#'
#' \strong{Label placement.}
#' \code{label.x.npc} / \code{label.y.npc} accept keywords or numeric positions
#' in \eqn{[0,1]}. Absolute positions via \code{label.x} / \code{label.y}
#' take precedence. The stat returns one row per panel; pass text aesthetics and
#' styling via \code{...}.
#'
#' \strong{Dependencies.}
#' Adjusted-mode center testing requires \pkg{quantreg}. If it is not available,
#' use the unadjusted mode (leave \code{confounder = NULL}) or ensure
#' \pkg{quantreg} is installed.
#'
#' @return
#' A \code{ggplot2} layer that draws a single label per panel. The underlying
#' computed data frame has columns:
#' \itemize{
#'   \item \code{x}, \code{y}: label coordinates
#'   \item \code{label}: combined multi-line p-value string
#'   \item \code{hjust}, \code{vjust}: justifications passed to the geom
#' }
#'
#' @examples
#' # Unadjusted example
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = factor(am), y = mpg)) +
#'   geom_boxplot() +
#'   stat_deshape_compare(size = 3.5)
#' p
#'
#' # Adjusted example (requires quantreg)
#' \donttest{
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   ggplot(mtcars, aes(x = factor(am), y = mpg)) +
#'     geom_boxplot() +
#'     stat_deshape_compare(
#'       confounder = c("hp", "wt"),
#'       confounder_data = mtcars,
#'       label.x.npc = "right",
#'       label.y.npc = "top",
#'       size = 3.5
#'     )
#' }
#' }
#'
#' @seealso
#' \code{\link{deshape_perm_multi}}, \code{\link{deshape_wald_contrast}},
#' \code{\link[ggplot2]{geom_text}}, \code{\link[ggplot2]{layer}}
#'
#' @importFrom ggplot2 layer aes ggproto Stat
#' @importFrom stats as.formula anova
#' @importFrom quantreg rq
#'
#' @export
stat_deshape_compare <- function(
    mapping = NULL,
    data = NULL,
    confounder = NULL,             # character vector of column names or NULL
    perm = 10000,
    seed = 2025,
    label.sep = "\n",              # separator between the three lines (default newline)
    label.prefix = c("Center", "Dispersion", "Asymmetry"),
    label.x.npc = "left",          # or numeric in [0,1]
    label.y.npc = "top",           # or numeric in [0,1]
    label.x = NULL,                # absolute data coordinates (overrides *_npc)
    label.y = NULL,                # absolute data coordinates (overrides *_npc)
    vjust = 0,
    geom = "text",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    confounder_data = NULL,
    ...
) {
  layer(
    stat = StatDeshapeCompare,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      confounder = confounder,
      perm = perm,
      seed = seed,
      label.sep = label.sep,
      label.prefix = label.prefix,
      label.x.npc = label.x.npc,
      label.y.npc = label.y.npc,
      label.x = label.x,
      label.y = label.y,
      vjust = vjust,
      na.rm = na.rm,
      confounder_data = confounder_data,
      ...
    )
  )
}

StatDeshapeCompare <- ggplot2::ggproto(
  "StatDeshapeCompare", ggplot2::Stat,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(hjust = ..hjust.., vjust = ..vjust..),
  
  compute_panel = function(data, scales,
                           confounder, perm, seed,
                           label.sep, label.prefix,
                           label.x.npc, label.y.npc,
                           label.x, label.y, vjust,
                           confounder_data = NULL) {
    
    # ---- helpers ------------------------------------------------------------
    fmt_p <- function(p) {
      # normalize to a single numeric if possible
      if (is.null(p) || length(p) == 0) return("NA")
      p <- p[1]
      if (!is.finite(p) || is.na(p)) return("NA")
      if (p < 2.2e-16) return("< 2.2e-16")
      as.character(signif(p, 4))
    }
    
    
    # Position helpers: support "left|right|center" / "top|bottom|center" or numeric [0,1]
    pos_from_npc <- function(npc, minv, maxv, axis = c("x", "y")) {
      axis <- match.arg(axis)
      if (is.null(npc)) return(list(v = NA_real_, just = 0.5))
      if (is.character(npc)) {
        npc <- tolower(npc)
        if (axis == "x") {
          if (npc %in% c("left"))   return(list(v = -Inf, just = -0.1))
          if (npc %in% c("right"))  return(list(v =  Inf, just =  1.1))
          if (npc %in% c("center","centre","middle")) return(list(v = (minv + maxv)/2, just = 0.5))
        } else {
          if (npc %in% c("top"))    return(list(v =  Inf, just =  1.1))
          if (npc %in% c("bottom")) return(list(v = -Inf, just = -0.1))
          if (npc %in% c("center","centre","middle")) return(list(v = (minv + maxv)/2, just = 0.5))
        }
        # fallback
        return(list(v = (minv + maxv)/2, just = 0.5))
      } else if (is.numeric(npc) && length(npc) == 1L && is.finite(npc)) {
        npc <- max(0, min(1, npc))
        return(list(v = minv + npc * (maxv - minv), just = 0.5))
      }
      list(v = (minv + maxv)/2, just = 0.5)
    }
    
    # Build working data.frame with generic column names to avoid relying on original names
    df <- data.frame(
      y = data$y,
      x = factor(data$x) # treat groups categorically
    )
    
    # If confounders were requested, try to carry them over from 'data'
    # If confounders were requested, pull them from confounder_data
    if (!is.null(confounder) && length(confounder) > 0) {
      if (is.null(confounder_data)) {
        stop("`confounder_data` must be supplied when `confounder` is not NULL.")
      }
      missing_cols <- setdiff(confounder, colnames(confounder_data))
      if (length(missing_cols) > 0) {
        stop("These confounder columns are not in `confounder_data`: ",
             paste(missing_cols, collapse = ", "))
      }
      
      # Align rows from the panel `data` to `confounder_data`.
      # Preferred: if rownames of `data` look like original row indices, use them.
      idx <- suppressWarnings(as.integer(rownames(data)))
      if (length(idx) == nrow(data) &&
          all(is.finite(idx)) &&
          all(idx >= 1) &&
          max(idx) <= nrow(confounder_data)) {
        
        df[confounder] <- confounder_data[idx, confounder, drop = FALSE]
        
      } else {
        # Fallback: assume `confounder_data` is in the same order as the data
        # used to build the layer and take the first n rows for this panel.
        # (For faceting or filtering, consider adding a stable row_id to align precisely.)
        if (nrow(confounder_data) < nrow(df)) {
          stop("`confounder_data` has fewer rows than the panel data; cannot align.")
        }
        df[confounder] <- confounder_data[seq_len(nrow(df)), confounder, drop = FALSE]
        warning("Row alignment fell back to sequential matching. ",
                "For faceting/filters, add a stable row id to align precisely.")
      }
    }
    
    
    # ---- compute p-values ---------------------------------------------------
    if (is.null(confounder) || length(confounder) == 0) {
      # Unadjusted: use your permutation engine for the three modes
      set.seed(seed)
      center_obj <- deshape_perm_multi(
        y ~ x, data = df, mode = "center", perm = perm, seed = seed
      )
      disp_obj <- deshape_perm_multi(
        y ~ x, data = df, mode = "dispersion", perm = perm, seed = seed
      )
      asym_obj <- deshape_perm_multi(
        y ~ x, data = df, mode = "skewness", perm = perm, seed = seed
      )
      
      p_center <- center_obj
      p_disp   <- disp_obj
      p_asym   <- asym_obj
      
    } else {
      # Adjusted: Center via median (tau=0.5) QR ANOVA; Disp/Asym via Wald-contrast
      rhs_full <- paste(c("x", confounder), collapse = " + ")
      f_full   <- stats::as.formula(paste("y ~", rhs_full))
      f_null   <- stats::as.formula(paste("y ~", paste(confounder, collapse = " + ")))
      
      # Center (location) via quantile regression ANOVA
      fit_full <- quantreg::rq(f_full, data = df, tau = 0.5, method = "fn")
      fit_null <- quantreg::rq(f_null, data = df, tau = 0.5, method = "fn")
      a <- stats::anova(fit_null, fit_full, joint = TRUE)
      # Extract p-value from the anova table (second row, last column typically)
      p_center <- tryCatch({
        # Try last column of second row
        as.numeric(a$table[1,4])
      }, error = function(e) NA_real_)
      
      # Dispersion
      disp_obj <- deshape_wald_contrast(
        f_full, data = df, mode = "dispersion",
        alternative = "two.sided", kernel = "gaussian"
      )
      p_disp <-disp_obj$p_value
      
      # Asymmetry (symmetry test)
      asym_obj <- deshape_wald_contrast(
        f_full, data = df, mode = "symmetry",
        alternative = "two.sided", kernel = "gaussian"
      )
      p_asym <- asym_obj$p_value
    }
    
    # ---- assemble label -----------------------------------------------------
    lab <- paste0(
      label.prefix[1], ": p = ", fmt_p(p_center), label.sep,
      label.prefix[2], ": p = ", fmt_p(p_disp),   label.sep,
      label.prefix[3], ": p = ", fmt_p(p_asym)
    )
    
    # ---- choose where to draw the label ------------------------------------
    xrange <- range(data$x, na.rm = TRUE)
    yrange <- range(data$y, na.rm = TRUE)
    
    if (!is.null(label.x)) {
      x_pos <- label.x
      h_just <- 0.5
    } else {
      x_calc <- pos_from_npc(label.x.npc, xrange[1], xrange[2], axis = "x")
      x_pos <- x_calc$v
      h_just <- x_calc$just
    }
    
    if (!is.null(label.y)) {
      y_pos <- label.y
      v_just <- 0.5
    } else {
      y_calc <- pos_from_npc(label.y.npc, yrange[1], yrange[2], axis = "y")
      y_pos <- y_calc$v
      v_just <- y_calc$just
    }
    
    # Return a one-row data frame; geom_text/label will use `label`
    data.frame(
      x = x_pos,
      y = y_pos,
      label = lab,
      hjust = h_just,
      vjust = v_just + vjust
    )
  }
)
