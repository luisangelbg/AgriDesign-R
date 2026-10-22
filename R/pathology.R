# Plant-pathology and field-plot tools: area under the disease progress curve
# (AUDPC) and stairs (AUDPS), absolute, relative and standardised, for one
# unit, many units or a long data frame of assessments; Smith's index of soil
# heterogeneity from a uniformity trial with the optimum plot size (Smith,
# Hatheway).

# ---------- AUDPC and AUDPS ----------
# AUDPC = sum (y_i + y_{i+1}) / 2 (t_{i+1} - t_i)         (trapezoids)
# AUDPS = AUDPC + (y_1 + y_n) / 2 * D / (n - 1)           (Simko and Piepho 2012)
# with D = t_n - t_1. The stairs extend half a mean interval before the first
# and after the last assessment, so the largest AUDPS is ymax D n / (n - 1):
# relative AUDPS divides by it, relative AUDPC by ymax D. The standardised
# areas divide by the width only (mean intensity, in the units of y).

.dpc_one <- function(t, y, ymax) {
  y <- as.numeric(y); t <- as.numeric(t)
  ok <- !is.na(y); t <- t[ok]; y <- y[ok]; n <- length(y)
  if (n < 2) return(c(n = n, duration = NA, audpc = NA, audps = NA, audpc_rel = NA, audps_rel = NA,
                      audpc_std = NA, audps_std = NA))
  D <- t[n] - t[1]
  a <- sum((y[-1] + y[-n]) / 2 * diff(t))
  s <- a + (y[1] + y[n]) / 2 * D / (n - 1)
  wS <- D * n / (n - 1)
  c(n = n, duration = D, audpc = a, audps = s, audpc_rel = a / (ymax * D), audps_rel = s / (ymax * wS),
    audpc_std = a / D, audps_std = s / wS)
}

.dpc_check_y <- function(y, ymax) {
  if (any(y < 0, na.rm = TRUE)) stop("disease intensities cannot be negative", call. = FALSE)
  if (length(ymax) != 1 || !is.finite(ymax) || ymax <= 0) stop("'ymax' must be a positive number", call. = FALSE)
  if (any(y > ymax, na.rm = TRUE))
    warning("some intensities exceed 'ymax' (", ymax, "): relative areas can be larger than 1", call. = FALSE)
}

.dpc_areas <- function(evaluation, dates, ymax) {
  if (is.data.frame(evaluation)) {
    if (!all(vapply(evaluation, is.numeric, NA))) stop("'evaluation' must be numeric", call. = FALSE)
    evaluation <- as.matrix(evaluation)
  }
  if (is.null(dim(evaluation))) evaluation <- matrix(evaluation, nrow = 1)
  if (!is.numeric(evaluation)) stop("'evaluation' must be numeric", call. = FALSE)
  if (inherits(dates, "Date")) dates <- as.numeric(dates)
  if (!is.numeric(dates) || anyNA(dates)) stop("'dates' must be numbers (days) or Date values without NA", call. = FALSE)
  if (length(dates) != ncol(evaluation))
    stop("give one date per assessment: ", ncol(evaluation), " assessments and ", length(dates), " dates", call. = FALSE)
  if (length(dates) < 2) stop("at least two assessments are needed", call. = FALSE)
  if (any(diff(dates) <= 0)) stop("'dates' must be strictly increasing", call. = FALSE)
  .dpc_check_y(evaluation, ymax)
  out <- t(apply(evaluation, 1, function(y) .dpc_one(dates, y, ymax)))
  if (any(out[, "n"] < 2)) warning("units with fewer than two assessments give NA", call. = FALSE)
  rownames(out) <- rownames(evaluation)
  out
}

#' Area under the disease progress curve and stairs
#'
#' `audpc()` integrates the disease progress curve with trapezoids (Shaner
#' and Finney 1977); `audps()` gives the area under the disease progress
#' stairs of Simko and Piepho (2012), which weighs the first and the last
#' assessment as much as the others: `AUDPS = AUDPC + (y1 + yn) / 2 *
#' D / (n - 1)`, with `D` the time from the first to the last assessment.
#'
#' `type = "relative"` divides the area by the largest possible one, so it
#' lies between 0 and 1: `ymax * D` for the AUDPC and `ymax * D * n / (n - 1)`
#' for the AUDPS, whose stairs are `n / (n - 1)` times as wide as the curve.
#' `type = "standardized"` divides by the width only and gives the mean
#' intensity over the epidemic, in the units of the assessments.
#' Missing assessments of a unit are skipped (its area then spans its own
#' assessments only).
#'
#' @param evaluation Disease intensities: a vector (one unit) or a matrix or
#'   data frame with one row per unit and one column per assessment.
#' @param dates Times of the assessments (days or [Date] values), strictly
#'   increasing.
#' @param type `"absolute"`, `"relative"` or `"standardized"`.
#' @param ymax Largest possible intensity (100 for percentages, 1 for
#'   proportions, the top of a rating scale).
#' @return A numeric vector with one area per unit.
#' @references Shaner, G. and Finney, R. E. (1977). The effect of nitrogen
#'   fertilization on the expression of slow-mildewing resistance in Knox
#'   wheat. *Phytopathology* 67, 1051-1056.
#'
#'   Simko, I. and Piepho, H.-P. (2012). The area under the disease progress
#'   stairs: calculation, advantage, and application. *Phytopathology* 102,
#'   381-389.
#' @seealso [audpc_by_unit()] to compute them for every plot of a trial.
#' @examples
#' days <- c(7, 14, 21, 28, 35)
#' sev <- c(1, 4, 12, 30, 52)
#' audpc(sev, days)
#' audps(sev, days)
#' audpc(sev, days, type = "relative")
#' audps(rbind(plot1 = sev, plot2 = sev / 2), days, type = "standardized")
#' @export
audpc <- function(evaluation, dates, type = c("absolute", "relative", "standardized"), ymax = 100) {
  type <- match.arg(type)
  a <- .dpc_areas(evaluation, dates, ymax)
  stats::setNames(unname(a[, switch(type, absolute = "audpc", relative = "audpc_rel", standardized = "audpc_std")]), rownames(a))
}

#' @rdname audpc
#' @export
audps <- function(evaluation, dates, type = c("absolute", "relative", "standardized"), ymax = 100) {
  type <- match.arg(type)
  a <- .dpc_areas(evaluation, dates, ymax)
  stats::setNames(unname(a[, switch(type, absolute = "audps", relative = "audps_rel", standardized = "audps_std")]), rownames(a))
}

#' AUDPC and AUDPS of every plot of a trial
#'
#' Takes the assessments in long format (one row per plot and date), sorts
#' them by time within each plot and returns one row per plot with its areas,
#' ready for [design_anova()].
#'
#' @param data A data frame in long format.
#' @param response Column with the disease intensity.
#' @param time Column with the time of assessment (numbers or [Date]).
#' @param units Column(s) that identify a plot, for instance
#'   `c("Block", "Treatment")`; they are kept in the result.
#' @inheritParams audpc
#' @return A data frame with the unit columns, the number of assessments
#'   (`n`), the time spanned (`duration`) and the areas `audpc`, `audps`,
#'   their relative (`audpc_rel`, `audps_rel`) and standardised (`audpc_std`,
#'   `audps_std`) forms.
#' @examples
#' a <- audpc_by_unit(blight_progress, "Severity_pct", "Day", c("Block", "Treatment"))
#' head(a)
#' fit <- design_anova(a, "audps", "Treatment", block = "Block")
#' compare_means(fit, "Treatment")
#' @export
audpc_by_unit <- function(data, response, time, units, ymax = 100) {
  if (!is.data.frame(data)) stop("'data' must be a data frame", call. = FALSE)
  .check_columns(data, c(response, time, units))
  d <- as.data.frame(data)[c(units, time, response)]
  if (!is.numeric(d[[response]])) stop("the response '", response, "' must be numeric", call. = FALSE)
  tm <- d[[time]]
  if (inherits(tm, "Date")) tm <- as.numeric(tm)
  if (!is.numeric(tm)) stop("the time column '", time, "' must hold numbers (days) or Date values", call. = FALSE)
  if (anyNA(tm) || anyNA(d[units])) stop("time and unit columns cannot have missing values", call. = FALSE)
  .dpc_check_y(d[[response]], ymax)
  key <- do.call(paste, c(lapply(d[units], as.character), sep = "\r"))
  ukeys <- unique(key)
  rows <- lapply(ukeys, function(u) {
    s <- which(key == u); o <- order(tm[s]); s <- s[o]
    if (any(duplicated(tm[s]))) stop("unit ", gsub("\r", " / ", u), " has two assessments at the same time", call. = FALSE)
    .dpc_one(tm[s], d[[response]][s], ymax)
  })
  out <- data.frame(d[match(ukeys, key), units, drop = FALSE], do.call(rbind, rows), row.names = NULL, check.names = FALSE)
  out$n <- as.integer(out$n)
  if (any(out$n < 2)) warning("units with fewer than two assessments give NA", call. = FALSE)
  if (length(unique(out$n)) > 1 || length(unique(out$duration)) > 1)
    warning("units differ in the number or span of assessments: compare the relative or standardised areas", call. = FALSE)
  out
}

# ---------- Smith's index of soil heterogeneity ----------

.divisors <- function(n) which(n %% seq_len(n) == 0)

#' Smith's index of soil heterogeneity
#'
#' From a uniformity trial (one crop, uniformly managed, harvested in small
#' basic units) builds every plot shape that tiles the field exactly
#' (`length` units along the rows times `width` units along the columns),
#' computes for each the variance per unit area among plots,
#' `V_x = var(T / x)` with `T` the plot totals and `x` the plot size in basic
#' units, and fits Smith's empirical law `V_x = V_1 / x^b` by least squares
#' on the log scale (weighted by the degrees of freedom of each variance,
#' `plots - 1`, unless `weights = "none"`).
#'
#' The index `b` lies usually between 0 and 1. Near 1, adjacent units are
#' nearly uncorrelated and larger plots reduce the variance per unit area
#' efficiently; near 0, adjacent units are much alike and enlarging plots
#' gains little, so more replicates pay off better. Given the ratio of the
#' cost per plot to the cost per unit area (`cost_ratio = K1 / K2`), Smith's
#' optimum plot size is `b K1 / ((1 - b) K2)` basic units; see
#' [hatheway_size()] for the size needed to detect a given difference.
#'
#' @param data A numeric matrix (or data frame) of yields of the basic
#'   units, with the rows and columns of the field.
#' @param cost_ratio Optional ratio of the cost per plot to the cost per unit
#'   area of experimental land (both per basic unit).
#' @param weights `"df"` (default) or `"none"`.
#' @param min_plots Smallest number of plots a shape must give to enter
#'   the fit (at least 2).
#' @param strips Use in the fit the shapes that span the whole length or
#'   width of the field. They are left out by default (as in agricolae):
#'   one strip per row or column averages out any gradient along it, so a
#'   few such shapes can dominate the fit, and they cannot be blocked.
#' @return An object of class `"ad_smith"` with the table of shapes
#'   (`table`: size, width, length, plots, df, variance, cv, used), the index
#'   `b`, the fitted `V1`, the R squared of the log-log fit, the mean yield,
#'   the coefficient of variation of the basic units (`cv1`) and Smith's
#'   optimum plot size (`optimum`).
#' @references Smith, H. F. (1938). An empirical law describing
#'   heterogeneity in the yields of agricultural crops. *Journal of
#'   Agricultural Science* 28, 1-23.
#'
#'   Gomez, K. A. and Gomez, A. A. (1984). *Statistical Procedures for
#'   Agricultural Research*, 2nd ed. Wiley, New York, chapter 12.
#' @examples
#' sm <- smith_index(uniformity_trial, cost_ratio = 2)
#' sm
#' plot(sm)
#' @export
smith_index <- function(data, cost_ratio = NULL, weights = c("df", "none"), min_plots = 2, strips = FALSE) {
  weights <- match.arg(weights)
  A <- if (is.data.frame(data)) {
    if (!all(vapply(data, is.numeric, NA))) stop("all columns of 'data' must be numeric yields", call. = FALSE)
    as.matrix(data)
  } else data
  if (!is.matrix(A) || !is.numeric(A)) stop("'data' must be a numeric matrix of basic units (field rows x columns)", call. = FALSE)
  if (anyNA(A)) stop("the uniformity trial cannot have missing units", call. = FALSE)
  if (min_plots < 2) stop("'min_plots' must be at least 2", call. = FALSE)
  nr <- nrow(A); nc <- ncol(A)
  if (nr * nc < 4) stop("the uniformity trial is too small", call. = FALSE)
  shapes <- expand.grid(length = .divisors(nr), width = .divisors(nc))
  shapes <- shapes[shapes$length * shapes$width < nr * nc, ]
  rows <- lapply(seq_len(nrow(shapes)), function(s) {
    l <- shapes$length[s]; w <- shapes$width[s]
    Tm <- rowsum(A, rep(seq_len(nr / l), each = l), reorder = FALSE)
    Tm <- t(rowsum(t(Tm), rep(seq_len(nc / w), each = w), reorder = FALSE))
    m <- as.vector(Tm) / (l * w)
    data.frame(size = l * w, width = w, length = l, plots = length(m), df = length(m) - 1,
               variance = stats::var(m), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$size, tab$width), ]
  rownames(tab) <- NULL
  mu <- mean(A)
  tab$cv <- sqrt(tab$variance) / mu * 100
  tab$used <- tab$plots >= min_plots & (strips | (tab$length < nr & tab$width < nc))
  u <- tab[tab$used, ]
  if (length(unique(u$size)) < 3) stop("too few plot sizes to fit Smith's law (need three or more)", call. = FALSE)
  if (any(u$variance <= 0)) stop("a plot shape has zero variance", call. = FALSE)
  w <- if (weights == "df") u$df else rep(1, nrow(u))
  X <- cbind(1, log(u$size)); ly <- log(u$variance)
  fw <- stats::lm.wfit(X, ly, w)
  cf <- unname(fw$coefficients)
  r2 <- 1 - sum(w * fw$residuals^2) / sum(w * (ly - stats::weighted.mean(ly, w))^2)
  b <- -cf[2]
  opt <- NULL
  if (!is.null(cost_ratio)) {
    if (length(cost_ratio) != 1 || !is.finite(cost_ratio) || cost_ratio <= 0)
      stop("'cost_ratio' must be a positive number", call. = FALSE)
    opt <- if (b > 0 && b < 1) b * cost_ratio / (1 - b) else NA_real_
  }
  structure(list(table = tab, b = b, V1 = exp(cf[1]), r_squared = r2, mean = mu, cv1 = tab$cv[tab$size == 1][1],
                 dim = c(rows = nr, cols = nc), weights = weights, cost_ratio = cost_ratio, optimum = opt,
                 min_plots = min_plots, strips = strips),
            class = "ad_smith")
}

#' @export
print.ad_smith <- function(x, digits = 4, ...) {
  cat("\nSmith's index of soil heterogeneity\n")
  cat(sprintf("Uniformity trial: %d x %d basic units, mean = %s, CV of the basic units = %.2f %%\n",
              x$dim[["rows"]], x$dim[["cols"]], format(x$mean, digits = digits), x$cv1))
  cat("\n"); print(x$table, row.names = FALSE, digits = digits)
  cat(sprintf("\nFit of V_x = V_1 / x^b on the log scale (%s, %d shapes with 'used' TRUE%s):\n  V_1 = %s, b = %.4f, R-squared = %.4f\n",
              if (x$weights == "df") "weighted by df" else "unweighted", sum(x$table$used),
              if (isTRUE(x$strips)) "" else "; strips spanning the whole field left out",
              format(x$V1, digits = digits), x$b, x$r_squared))
  cat(if (x$b >= 0.7) "b near 1: adjacent units are little correlated; larger plots reduce the variance per unit area efficiently.\n"
      else if (x$b <= 0.3) "b near 0: adjacent units are much alike; enlarging plots gains little, more replicates pay off better.\n"
      else "Intermediate b: adjacent units are moderately correlated.\n")
  if (!is.null(x$cost_ratio)) {
    if (is.na(x$optimum)) cat("Smith's optimum plot size is undefined when b is not between 0 and 1.\n")
    else cat(sprintf("Optimum plot size (Smith, cost per plot / cost per unit area = %g): %s basic units\n",
                     x$cost_ratio, format(x$optimum, digits = digits)))
  }
  invisible(x)
}

#' @export
plot.ad_smith <- function(x, col = "#2f7d4f", main = "Smith's law of soil heterogeneity",
                          xlab = "Plot size (basic units)", ylab = "Variance per unit area", ...) {
  tb <- x$table
  op <- .style_par(); on.exit(graphics::par(op))
  xs <- exp(seq(log(min(tb$size)), log(max(tb$size)), length.out = 100))
  ys <- x$V1 / xs^x$b
  graphics::plot(tb$size, tb$variance, log = "xy", type = "n", xlab = xlab, ylab = ylab,
                 ylim = range(c(tb$variance, ys)), ...)
  graphics::abline(h = pretty(range(tb$variance)), col = "#ececec", lwd = 0.8)
  graphics::lines(xs, ys, col = col, lwd = 2.2)
  graphics::points(tb$size, tb$variance, pch = 21, bg = ifelse(tb$used, col, "#ffffff"), col = .darken(col),
                   cex = 0.8 + 1.2 * sqrt(tb$df / max(tb$df)))
  graphics::legend("topright", bty = "n", cex = 0.85,
                   legend = sprintf("V = %s / x^%.3f   R%s = %.3f", format(signif(x$V1, 4)), x$b, .sup[1], x$r_squared))
  graphics::title(main = main, font.main = 1, cex.main = 1.05)
  .note(paste0("b = ", formatC(x$b, format = "f", digits = 3), .sep, "point size: df", .sep, "open: not in the fit",
               if (!is.null(x$optimum) && !is.na(x$optimum)) paste0(.sep, "optimum = ", format(x$optimum, digits = 3), " units") else ""))
  invisible(tb)
}

#' Plot size needed to detect a difference (Hatheway)
#'
#' Hatheway's (1961) relation between plot size, replication and the
#' difference that can be detected, from Smith's index:
#' `x^b = 2 (t1 + t2)^2 CV1^2 / (r d^2)`, where `CV1` is the coefficient of
#' variation of the basic units (%), `d` the true difference to detect as a
#' percentage of the mean, `r` the number of replicates, `t1` the two-sided
#' t for `alpha` and `t2` the one-sided t for the power wanted.
#'
#' @param index An `"ad_smith"` object from [smith_index()], or the value of
#'   `b` (then give `cv`).
#' @param reps Number(s) of replicates.
#' @param difference Difference(s) to detect, as a percentage of the mean.
#' @param alpha Significance level (two-sided).
#' @param power Probability of detecting the difference.
#' @param treatments Optional number of treatments: the t values then use
#'   the error degrees of freedom of a randomised complete block design,
#'   `(reps - 1) (treatments - 1)`; otherwise normal quantiles are used.
#' @param cv Coefficient of variation (%) of the basic units; by default
#'   the one of the uniformity trial.
#' @return A data frame with `reps`, `difference`, the error `df` and the
#'   plot size in basic units (`plot_size`).
#' @references Hatheway, W. H. (1961). Convenient plot size. *Agronomy
#'   Journal* 53, 279-280.
#' @examples
#' sm <- smith_index(uniformity_trial)
#' hatheway_size(sm, reps = 3:6, difference = c(10, 15), treatments = 8)
#' @export
hatheway_size <- function(index, reps, difference, alpha = 0.05, power = 0.80, treatments = NULL, cv = NULL) {
  if (inherits(index, "ad_smith")) { b <- index$b; if (is.null(cv)) cv <- index$cv1 }
  else b <- index
  if (length(b) != 1 || !is.finite(b) || b <= 0) stop("Smith's index 'b' must be a positive number", call. = FALSE)
  if (is.null(cv) || length(cv) != 1 || !(cv > 0)) stop("give the coefficient of variation 'cv' (%) of the basic units", call. = FALSE)
  if (any(reps < 2)) stop("'reps' must be at least 2", call. = FALSE)
  if (any(difference <= 0)) stop("'difference' must be positive (percentage of the mean)", call. = FALSE)
  .check_alpha(alpha); if (!(power > 0 && power < 1)) stop("'power' must be between 0 and 1", call. = FALSE)
  out <- expand.grid(reps = reps, difference = difference)
  out$df <- if (is.null(treatments)) Inf else (out$reps - 1) * (treatments - 1)
  t1 <- stats::qt(1 - alpha / 2, out$df); t2 <- stats::qt(power, out$df)
  out$plot_size <- (2 * (t1 + t2)^2 * cv^2 / (out$reps * out$difference^2))^(1 / b)
  out
}

# ---------- example data ----------

#' Late blight progress under four fungicide programmes (simulated)
#'
#' Simulated severity of late blight (percentage of leaf area affected) in
#' a potato trial with four fungicide programmes in four blocks, assessed
#' weekly from day 21 to day 63 after planting. The values were generated
#' with logistic epidemics (data-raw/make-data-pathology.R).
#'
#' @format A data frame with 112 rows and 4 columns:
#' \describe{
#'   \item{Block}{Block (1-4).}
#'   \item{Treatment}{Fungicide programme: Untreated, Contact, Systemic, Alternated.}
#'   \item{Day}{Days after planting.}
#'   \item{Severity_pct}{Disease severity, percent.}
#' }
#' @examples
#' a <- audpc_by_unit(blight_progress, "Severity_pct", "Day", c("Block", "Treatment"))
#' design_anova(a, "audpc", "Treatment", block = "Block")
"blight_progress"

#' Uniformity trial of 24 x 16 basic units (simulated)
#'
#' Simulated yields (kg) of a uniformity trial harvested in 384 basic units
#' (24 rows by 16 columns), with a fertility gradient along the rows and
#' patchy variation; suited to Smith's index of soil heterogeneity.
#'
#' @format A numeric matrix with 24 rows (`R1`...`R24`) and 16 columns
#'   (`C1`...`C16`).
#' @examples
#' smith_index(uniformity_trial)
"uniformity_trial"

#' Sorghum hybrids in a balanced incomplete block design (simulated)
#'
#' Simulated grain yields of seven sorghum hybrids in seven blocks of three
#' plots: every hybrid appears three times and every pair of hybrids meets
#' once in a block (t = 7, k = 3, r = 3, lambda = 1).
#'
#' @format A data frame with 21 rows and 3 columns:
#' \describe{
#'   \item{Block}{Block (1-7).}
#'   \item{Hybrid}{Hybrid (H1-H7).}
#'   \item{Yield_t_ha}{Grain yield, t/ha.}
#' }
#' @examples
#' durbin_test(bib_sorghum, "Yield_t_ha", "Hybrid", "Block")
"bib_sorghum"
