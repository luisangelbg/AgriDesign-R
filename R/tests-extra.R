# Further mean-separation and rank-based tests: the Waller-Duncan k-ratio t
# test (Bayesian), Durbin's test for balanced incomplete blocks, the van der
# Waerden normal-scores test and Mood's median test, each with its pairwise
# comparisons and compact letters.

# ---------- Waller-Duncan critical t ----------
# Model (Waller and Duncan 1969): treatment means ybar_i ~ N(mu_i, s2 / n),
# mu_i ~ N(mu, s2_mu) exchangeable; error MS with f df and treatment MS with q
# df. Given sigma^2 and sigma_T^2 = sigma^2 + n s2_mu, the posterior of a
# difference delta = mu_i - mu_j is N(beta d, 2 beta sigma^2 / n) with the
# shrinkage beta = 1 - sigma^2 / sigma_T^2. With the linear losses k0 delta
# (missing a positive difference) and k1 |delta| (declaring a wrong sign), the
# Bayes rule declares delta > 0 when E[delta+] / E[delta-] > K = k1 / k0.
#
# Writing t = d / sqrt(2 MSE / n), F the observed ratio and X ~ F(q, f) the
# ratio of the two chi-square variates (restricted to beta = 1 - X / F > 0),
# the error variance integrates out in closed form: given X, the scaled
# precision (f + q X) MSE / sigma^2 is chi-square with nu = f + q df, which
# turns the normal partial expectations into Student t ones:
#   E[delta+ | X] ~ t beta P_nu(theta sqrt(nu)) + sqrt(beta (f + q X)) Q_nu(theta)
#   E[delta- | X] ~ -t beta P_nu(-theta sqrt(nu)) + sqrt(beta (f + q X)) Q_nu(theta)
# with theta = t sqrt(beta / (f + q X)), P_nu the t distribution function and
# Q_nu(theta) = Gamma((nu - 1) / 2) / (2 sqrt(pi) Gamma(nu / 2)) (1 + theta^2)^(-(nu - 1) / 2).
# Both are then integrated over X in (0, F) with the F(q, f) density, and the
# critical t is the root of log(num / den) = log(K).

.wd_parts <- function(t, q, f, Fv) {
  nu <- f + q
  lc <- lgamma((nu - 1) / 2) - lgamma(nu / 2) - log(2 * sqrt(pi))
  # integrated on the probability scale P = pf(X, q, f): the F density is
  # absorbed and the integrand stays smooth for any q and any observed F
  inner <- function(P, sgn) {
    X <- stats::qf(P, q, f)
    beta <- pmax(0, 1 - X / Fv); h <- f + q * X
    th <- t * sqrt(beta / h)
    Q <- exp(lc - (nu - 1) / 2 * log1p(th^2))
    v <- sqrt(beta * h) * Q + sgn * t * beta * stats::pt(sgn * th * sqrt(nu), nu)
    v[!is.finite(v)] <- 0
    v
  }
  upper <- stats::pf(Fv, q, f)
  num <- stats::integrate(inner, 0, upper, sgn = 1, rel.tol = 1e-10, subdivisions = 1000L, stop.on.error = FALSE)$value
  den <- stats::integrate(inner, 0, upper, sgn = -1, rel.tol = 1e-10, subdivisions = 1000L, stop.on.error = FALSE)$value
  c(num = num, den = den)
}

.wd_crit <- function(K, q, f, Fv) {
  f <- min(f, 1e7)
  h <- function(t) {
    p <- .wd_parts(t, q, f, Fv)
    if (!is.finite(p[["den"]]) || p[["den"]] <= 0) return(50)
    log(p[["num"]]) - log(p[["den"]]) - log(K)
  }
  hi <- 5
  while (h(hi) < 0 && hi < 500) hi <- hi * 2
  stats::uniroot(h, c(1e-3, hi), tol = 1e-10)$root
}

#' Critical t of the Waller-Duncan k-ratio test
#'
#' The Bayesian critical value t_B(K, F, q, f) of Waller and Duncan (1969),
#' computed by numerical integration of the posterior expected losses (see
#' [waller_duncan()]). Two means differ when their difference exceeds
#' `t_B * SE`, the Bayes least significant difference. The value decreases
#' as the observed F grows: with a small F the means are shrunk towards
#' each other and few differences are declared.
#'
#' @param K Ratio of the seriousness of a type I error to that of a type II
#'   error (50, 100 and 500 roughly match alpha = 0.10, 0.05 and 0.01).
#' @param q,f Degrees of freedom of the treatments and of the error.
#' @param F Observed F ratio of the treatments.
#' @return A numeric vector of critical values (arguments are recycled).
#' @references Waller, R. A. and Duncan, D. B. (1969). A Bayes rule for the
#'   symmetric multiple comparisons problem. *Journal of the American
#'   Statistical Association* 64, 1484-1503.
#'
#'   Waller, R. A. and Kemp, K. E. (1975). Computations of Bayesian t-values
#'   for multiple comparisons. *Journal of Statistical Computation and
#'   Simulation* 4, 169-171.
#' @examples
#' waller_critical(100, q = 5, f = 24, F = 14.37)
#' waller_critical(100, q = 4, f = 20, F = c(2, 5, 10, 100))
#' @export
waller_critical <- function(K = 100, q, f, F) {
  if (any(!is.finite(K)) || any(K <= 1)) stop("'K' must be greater than 1", call. = FALSE)
  if (any(q < 1) || any(f < 1)) stop("degrees of freedom 'q' and 'f' must be at least 1", call. = FALSE)
  if (any(!(F > 0))) stop("'F' must be positive", call. = FALSE)
  unname(mapply(.wd_crit, K, q, f, F))
}

#' Waller-Duncan k-ratio t test
#'
#' Bayesian separation of the means of a treatment factor (Waller and Duncan
#' 1969). Instead of a fixed significance level, the user gives the ratio
#' `K` of the seriousness of a type I to a type II error; the critical t
#' then adapts to the evidence of the analysis of variance: with a large F
#' it approaches an unprotected LSD, with an F near 1 almost no difference is
#' declared. The result is an `"ad_comparison"` object, so [print()],
#' [plot()] and [plot_differences()] work as with [compare_means()].
#'
#' The F, its degrees of freedom and the error come from the stratum that
#' tests the factor in `fit` (for instance Error a for the main plots of a
#' split-plot). With unequal replication every pair uses its own standard
#' error, as in the Tukey-Kramer form.
#'
#' @inheritParams compare_means
#' @param factor One treatment factor of `fit`.
#' @param K Seriousness ratio of type I to type II errors (default 100).
#' @param alpha Optional nominal level that sets `K` by the usual
#'   correspondence: 0.10 gives K = 50, 0.05 gives K = 100 and 0.01 gives
#'   K = 500. When `NULL`, the `alpha` reported is the comparison-wise error
#'   rate implied by the critical t.
#' @return An object of class `"ad_comparison"` (see [compare_means()]),
#'   with the extra elements `K`, `t_critical` and `F`.
#' @references Waller, R. A. and Duncan, D. B. (1969). A Bayes rule for the
#'   symmetric multiple comparisons problem. *Journal of the American
#'   Statistical Association* 64, 1484-1503.
#' @seealso [waller_critical()], [compare_means()]
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' wd <- waller_duncan(fit, "Nitrogen_kg_ha")
#' wd
#' plot(wd)
#' @export
waller_duncan <- function(fit, factor, K = 100, alpha = NULL) {
  .check_fit(fit)
  if (length(factor) != 1) stop("the Waller-Duncan test compares the levels of one factor", call. = FALSE)
  .check_factors(fit, factor)
  if (!is.null(alpha)) {
    K <- switch(as.character(alpha), "0.1" = 50, "0.05" = 100, "0.01" = 500,
                stop("'alpha' must be 0.10, 0.05 or 0.01 (K = 50, 100, 500); otherwise give 'K'", call. = FALSE))
  }
  if (length(K) != 1 || !is.finite(K) || K <= 1) stop("'K' must be a single number greater than 1", call. = FALSE)
  row <- fit$table[fit$table$term == factor, ]
  if (!nrow(row) || is.na(row$F)) stop("the analysis of variance has no F test for '", factor, "'", call. = FALSE)
  lsm <- .lsm(fit, factor)
  residual <- lsm$error == "Residuals"
  err <- if (residual) list(ms = fit$residual[["ms"]], df = fit$residual[["df"]], label = fit$strata$Residuals$label)
         else fit$strata[[lsm$error]]
  tb <- waller_critical(K, row$df, row$df_error, row$F)
  g <- lsm$means; k <- nrow(g)
  se_diff <- .se_diff_fun(fit, lsm, err$ms, residual)
  ij <- .pair_index(k)
  diff <- g$mean[ij[, 1]] - g$mean[ij[, 2]]
  se <- vapply(seq_len(nrow(ij)), function(q) se_diff(ij[q, 1], ij[q, 2]), 0)
  tval <- diff / se
  pairs <- data.frame(i = ij[, 1], j = ij[, 2], level_i = g$label[ij[, 1]], level_j = g$label[ij[, 2]],
                      diff = diff, se = se, statistic = tval, p = NA_real_, p_adj = NA_real_,
                      different = abs(tval) > tb, critical = tb * se, stringsAsFactors = FALSE)
  bal <- length(unique(g$n)) == 1
  g$group <- .letters(g$mean, .diff_matrix(k, ij, ifelse(pairs$different, 0, 1), 0.5))
  if (!is.null(fit$transform)) g$back_transformed <- fit$transform$inv(g$mean)
  a_cw <- 2 * stats::pt(-tb, err$df)
  note <- sprintf("Bayesian k-ratio rule: K = %g, F = %.4g on %g and %g df, critical t = %.4f (comparison-wise alpha about %.3f).%s",
                  K, row$F, row$df, row$df_error, tb, a_cw,
                  if (bal) "" else " Unequal replication: own SE for each pair.")
  structure(list(means = g, pairs = pairs, msd = if (bal) c(value = tb * se_diff(1, 2)) else NULL,
                 note = note, method = "waller", method_name = sprintf("Waller-Duncan k-ratio t test (K = %g)", K),
                 alpha = if (is.null(alpha)) signif(a_cw, 3) else alpha,
                 error = list(ms = err$ms, df = err$df, label = err$label),
                 response = fit$response, transform = fit$transform, factor = factor, by = NULL,
                 obs = data.frame(label = as.character(fit$data[[factor]]), y = fit$y, stringsAsFactors = FALSE),
                 K = K, t_critical = tb, F = row$F),
            class = "ad_comparison")
}

# ---------- shared input of the rank-based tests ----------

# y, treatment factor g and (optionally) block factor from a fit or a data frame.
.rank_input <- function(x, response, treatment, block, what, blocks = c("none", "required", "optional")) {
  blocks <- match.arg(blocks)
  if (inherits(x, "ad_anova")) {
    r <- x$roles
    if (length(r$treatment) != 1) stop(what, " compares the levels of one treatment factor", call. = FALSE)
    if (!is.null(r$covariate)) stop(what, " does not handle covariates", call. = FALSE)
    if (!is.null(r$row) || !is.null(r$col)) stop(what, " does not handle row-column designs", call. = FALSE)
    y <- x$y; g <- x$data[[r$treatment]]
    blk <- if (length(r$block)) interaction(x$data[r$block], drop = TRUE, sep = ":") else NULL
    response <- x$response; treatment <- r$treatment
  } else if (is.data.frame(x)) {
    if (is.null(response) || is.null(treatment))
      stop("with a data frame give 'response' and 'treatment' (column names)", call. = FALSE)
    if (length(response) != 1 || length(treatment) != 1) stop("give one 'response' and one 'treatment' column", call. = FALSE)
    .check_columns(x, c(response, treatment, block))
    d <- as.data.frame(x)[c(response, treatment, block)]
    if (!is.numeric(d[[response]])) stop("the response '", response, "' must be numeric", call. = FALSE)
    d <- d[stats::complete.cases(d), , drop = FALSE]
    y <- d[[response]]
    g <- factor(as.character(d[[treatment]]), levels = .levels_of(d[[treatment]]))
    blk <- if (length(block)) interaction(lapply(d[block], as.character), drop = TRUE, sep = ":") else NULL
  } else stop("'x' must be a data frame or come from design_anova()", call. = FALSE)
  if (blocks == "none" && !is.null(blk))
    stop(what, " is for a one-way layout without blocks; use durbin_test() or nonparametric() for blocked designs", call. = FALSE)
  if (blocks == "required" && is.null(blk)) stop(what, " needs a block column", call. = FALSE)
  g <- droplevels(g)
  if (nlevels(g) < 2) stop("the treatment needs at least two levels", call. = FALSE)
  list(y = y, g = g, blk = blk, response = response, factor = treatment)
}

.rank_pairs <- function(lv, ij, stat_name, stat, p, p_adjust, alpha) {
  padj <- stats::p.adjust(p, p_adjust)
  data.frame(level_i = lv[ij[, 1]], level_j = lv[ij[, 2]], statistic_name = stat_name,
             statistic = unname(stat), p = unname(p), p_adj = unname(padj), different = unname(padj < alpha),
             stringsAsFactors = FALSE)
}

.rank_result <- function(cls, method, test, table, pairs, groups, center, posthoc, p_adjust, alpha, inp, extra) {
  k <- nrow(groups)
  groups$group <- .letters(center, .diff_matrix(k, .pair_index(k), pairs$p_adj, alpha))
  groups <- groups[order(-center), ]
  rownames(groups) <- NULL
  structure(c(list(method = method, test = test, table = table, pairs = pairs, groups = groups,
                   posthoc = posthoc, p_adjust = p_adjust, alpha = alpha, factor = inp$factor,
                   response = inp$response, available = method), extra),
            class = c(cls, "ad_nonparametric"))
}

.check_alpha <- function(alpha) {
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a number between 0 and 1", call. = FALSE)
}

# ---------- Durbin ----------

#' Durbin test for balanced incomplete blocks
#'
#' Rank-based comparison of treatments in a balanced incomplete block design
#' (every block has `k` of the `t` treatments, every treatment appears `r`
#' times): the observations are ranked within each block and the rank sums
#' are compared (Durbin 1951). With complete blocks it is Friedman's test.
#' Pairwise comparisons use Conover's t test on the rank sums with
#' `b k - b - t + 1` degrees of freedom; ties get average ranks.
#'
#' @param x An `"ad_anova"` object from [design_anova()] with one treatment
#'   factor and blocks, or a data frame.
#' @param response,treatment,block Column names when `x` is a data frame;
#'   several block columns (replicate, block within replicate) are combined.
#' @param alpha Significance level of the pairwise comparisons.
#' @param p_adjust Adjustment of the pairwise p-values, as in
#'   [stats::p.adjust()]; `"none"` gives Conover's protected comparisons
#'   (to be read only when the overall test is significant).
#' @return An object of class `c("ad_durbin", "ad_nonparametric")` with the
#'   test (chi-square and its F form), the pairwise comparisons, the groups
#'   with their rank sums and letters, the design parameters (`design`) and
#'   the least significant difference of rank sums (`lsd_ranks`).
#' @references Durbin, J. (1951). Incomplete blocks in ranking experiments.
#'   *British Journal of Statistical Psychology* 4, 85-90.
#'
#'   Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd ed.
#'   Wiley, New York, section 5.9.
#' @examples
#' durbin_test(bib_sorghum, "Yield_t_ha", "Hybrid", "Block")
#' @export
durbin_test <- function(x, response = NULL, treatment = NULL, block = NULL, alpha = 0.05, p_adjust = "holm") {
  .check_alpha(alpha)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  inp <- .rank_input(x, response, treatment, block, "Durbin's test", "required")
  y <- inp$y; g <- inp$g; blk <- droplevels(inp$blk)
  tab <- table(blk, g)
  if (any(tab > 1)) stop("a treatment appears more than once in a block; Durbin's test needs one plot per block and treatment", call. = FALSE)
  ks <- rowSums(tab); rs <- colSums(tab)
  if (length(unique(ks)) != 1 || length(unique(rs)) != 1)
    stop("Durbin's test needs blocks of equal size and treatments with equal replication (a balanced incomplete block design)", call. = FALSE)
  ntr <- ncol(tab); b <- nrow(tab); k <- ks[[1]]; r <- rs[[1]]
  if (k < 2) stop("blocks must hold at least two plots", call. = FALSE)
  conc <- crossprod(tab)[upper.tri(diag(ntr))]
  if (length(unique(conc)) != 1)
    warning("pairs of treatments do not meet equally often in the blocks: the design is not balanced and the test is approximate", call. = FALSE)
  rk <- stats::ave(y, blk, FUN = rank)
  R <- as.numeric(tapply(rk, g, sum))
  A <- sum(rk^2); C <- b * k * (k + 1)^2 / 4
  if (A - C <= 0) stop("all observations are tied within every block", call. = FALSE)
  T1 <- (ntr - 1) * sum((R - r * (k + 1) / 2)^2) / (A - C)
  df2 <- b * k - b - ntr + 1
  if (df2 < 1) stop("too few blocks for Durbin's test", call. = FALSE)
  T2 <- (T1 / (ntr - 1)) / ((b * (k - 1) - T1) / df2)
  se <- sqrt(max(0, (A - C) * 2 * r / df2 * (1 - T1 / (b * (k - 1)))))
  ij <- .pair_index(ntr); lv <- levels(g)
  tt <- (R[ij[, 1]] - R[ij[, 2]]) / se
  pairs <- .rank_pairs(lv, ij, "t", tt, 2 * stats::pt(-abs(tt), df2), p_adjust, alpha)
  groups <- data.frame(level = lv, n = as.integer(tabulate(g, ntr)), median = as.numeric(tapply(y, g, stats::median)),
                       mean = as.numeric(tapply(y, g, mean)), mean_rank = R / r, rank_sum = R, stringsAsFactors = FALSE)
  test <- list(name = "Durbin rank test for balanced incomplete blocks", statistic = c(chisq = T1), df = ntr - 1,
               p = stats::pchisq(T1, ntr - 1, lower.tail = FALSE), effect = NULL,
               F_form = c(F = T2, df1 = ntr - 1, df2 = df2,
                          p = if (is.finite(T2) && T2 >= 0) stats::pf(T2, ntr - 1, df2, lower.tail = FALSE) else 0))
  .rank_result("ad_durbin", "durbin", test, NULL, pairs, groups, R, "Conover t on rank sums", p_adjust, alpha, inp,
               list(design = c(treatments = ntr, blocks = b, k = k, r = r, lambda = r * (k - 1) / (ntr - 1)),
                    lsd_ranks = stats::qt(1 - alpha / 2, df2) * se, df_error = df2))
}

#' @export
print.ad_durbin <- function(x, digits = 4, ...) {
  print.ad_nonparametric(x, digits = digits)
  d <- x$design; f <- x$test$F_form
  cat(sprintf("\nDesign: %g treatments in %g blocks of %g plots, r = %g, lambda = %s\n", d[["treatments"]], d[["blocks"]],
              d[["k"]], d[["r"]], format(d[["lambda"]], digits = digits)))
  cat(sprintf("F form: F = %s on %g and %g df, p = %s\n", format(f[["F"]], digits = digits), f[["df1"]], f[["df2"]],
              format.pval(f[["p"]], digits = digits)))
  if (x$p_adjust == "none")
    cat(sprintf("Least significant difference of rank sums = %s (alpha = %g)\n", format(x$lsd_ranks, digits = digits), x$alpha))
  invisible(x)
}

# ---------- van der Waerden ----------

#' van der Waerden normal-scores test
#'
#' Compares the treatments of a one-way layout (no blocks) after replacing
#' every observation by its normal score `qnorm(R / (N + 1))`, where `R` is
#' its rank among all `N` observations (average ranks for ties). The
#' statistic is chi-square with `t - 1` degrees of freedom; the pairwise
#' comparisons are Conover's t tests on the mean scores with `N - t`
#' degrees of freedom. It is as powerful as the F test under normality and
#' more robust than it.
#'
#' @inheritParams durbin_test
#' @param x An `"ad_anova"` object of a one-factor design without blocks,
#'   or a data frame.
#' @param response,treatment Column names when `x` is a data frame.
#' @return An object of class `c("ad_waerden", "ad_nonparametric")`; the
#'   groups carry the mean normal score of every level (`mean_score`).
#' @references van der Waerden, B. L. (1952). Order tests for the
#'   two-sample problem and their power. *Indagationes Mathematicae* 14,
#'   453-458.
#'
#'   Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd ed.
#'   Wiley, New York, section 5.10.
#' @examples
#' waerden_test(bean_varieties, "Yield_kg_plot", "Variety")
#' @export
waerden_test <- function(x, response = NULL, treatment = NULL, alpha = 0.05, p_adjust = "holm") {
  .check_alpha(alpha)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  inp <- .rank_input(x, response, treatment, NULL, "The van der Waerden test", "none")
  y <- inp$y; g <- inp$g
  N <- length(y); k <- nlevels(g); n <- tabulate(g, k)
  if (N - k < 1) stop("too few observations: every level needs replicates", call. = FALSE)
  A <- stats::qnorm(rank(y) / (N + 1))
  s2 <- sum(A^2) / (N - 1)
  if (s2 <= 0) stop("all observations are tied", call. = FALSE)
  Ab <- as.numeric(tapply(A, g, mean))
  T1 <- sum(n * Ab^2) / s2
  ij <- .pair_index(k); lv <- levels(g)
  se <- sqrt(s2 * (N - 1 - T1) / (N - k) * (1 / n[ij[, 1]] + 1 / n[ij[, 2]]))
  tt <- (Ab[ij[, 1]] - Ab[ij[, 2]]) / se
  pairs <- .rank_pairs(lv, ij, "t", tt, 2 * stats::pt(-abs(tt), N - k), p_adjust, alpha)
  groups <- data.frame(level = lv, n = n, median = as.numeric(tapply(y, g, stats::median)),
                       mean = as.numeric(tapply(y, g, mean)), mean_rank = as.numeric(tapply(rank(y), g, mean)),
                       mean_score = Ab, stringsAsFactors = FALSE)
  test <- list(name = "van der Waerden normal-scores test", statistic = c(chisq = T1), df = k - 1,
               p = stats::pchisq(T1, k - 1, lower.tail = FALSE), effect = NULL)
  .rank_result("ad_waerden", "waerden", test, NULL, pairs, groups, Ab, "Conover t on mean normal scores", p_adjust,
               alpha, inp, list(score_variance = s2, df_error = N - k))
}

#' @export
print.ad_waerden <- function(x, digits = 4, ...) {
  print.ad_nonparametric(x, digits = digits)
  cat(sprintf("\nVariance of the normal scores = %s; pairwise t with %g df\n", format(x$score_variance, digits = digits), x$df_error))
  invisible(x)
}

# ---------- Mood's median test ----------

.chisq_stat <- function(O, correct) {
  E <- outer(rowSums(O), colSums(O)) / sum(O)
  if (any(E == 0)) return(0)
  d <- abs(O - E)
  if (correct && all(dim(O) == 2)) d <- d - pmin(0.5, d)
  sum(d^2 / E)
}

#' Mood's median test
#'
#' Compares the treatments of a one-way layout by the number of observations
#' above the grand median: a chi-square test on the 2 x t table of counts
#' above and at-or-below the median (Mood 1950; Conover 1999). Pairs of
#' treatments are compared with the same test on their 2 x 2 table, either
#' around the median of the two treatments together (`pair_median = "pair"`,
#' the default) or around the grand median (`"grand"`), or with Fisher's
#' exact test (`exact = TRUE`). It assumes little more than independence,
#' and resists outliers and unequal spreads, at the cost of power.
#'
#' @inheritParams waerden_test
#' @param correct Yates' continuity correction in the 2 x 2 tables (the
#'   overall test is corrected only when there are two treatments).
#' @param pair_median Median used for the pairwise tables.
#' @param exact Use Fisher's exact test for the pairwise tables.
#' @return An object of class `c("ad_median_test", "ad_nonparametric")`;
#'   `table` holds the counts above and at-or-below the grand median and
#'   `grand_median` the median of all observations.
#' @references Mood, A. M. (1950). *Introduction to the Theory of
#'   Statistics*. McGraw-Hill, New York, 394-399.
#'
#'   Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd ed.
#'   Wiley, New York, section 4.3.
#' @examples
#' median_test(bean_varieties, "Yield_kg_plot", "Variety")
#' @export
median_test <- function(x, response = NULL, treatment = NULL, alpha = 0.05, p_adjust = "holm",
                        correct = FALSE, pair_median = c("pair", "grand"), exact = FALSE) {
  .check_alpha(alpha)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  pair_median <- match.arg(pair_median)
  inp <- .rank_input(x, response, treatment, NULL, "The median test", "none")
  y <- inp$y; g <- inp$g
  k <- nlevels(g); lv <- levels(g)
  M <- stats::median(y)
  above <- as.numeric(tapply(y > M, g, sum)); n <- tabulate(g, k)
  O <- rbind(above, n - above)
  X2 <- .chisq_stat(O, correct && k == 2)
  ij <- .pair_index(k)
  pr <- t(vapply(seq_len(nrow(ij)), function(q) {
    sel <- g %in% lv[ij[q, ]]
    yy <- y[sel]; gg <- factor(g[sel], levels = lv[ij[q, ]])
    Mp <- if (pair_median == "pair") stats::median(yy) else M
    a <- as.numeric(tapply(yy > Mp, gg, sum)); nn <- tabulate(gg, 2)
    O2 <- rbind(a, nn - a)
    if (exact) c(NA_real_, stats::fisher.test(O2)$p.value, Mp)
    else { s <- .chisq_stat(O2, correct); c(s, stats::pchisq(s, 1, lower.tail = FALSE), Mp) }
  }, c(0, 0, 0)))
  pairs <- .rank_pairs(lv, ij, if (exact) "Fisher" else "chisq", pr[, 1], pr[, 2], p_adjust, alpha)
  pairs$median <- pr[, 3]
  groups <- data.frame(level = lv, n = n, median = as.numeric(tapply(y, g, stats::median)),
                       mean = as.numeric(tapply(y, g, mean)), mean_rank = as.numeric(tapply(rank(y), g, mean)),
                       above = above, stringsAsFactors = FALSE)
  test <- list(name = "Mood's median test", statistic = c(chisq = X2), df = k - 1,
               p = stats::pchisq(X2, k - 1, lower.tail = FALSE), effect = NULL)
  tab <- data.frame(level = lv, n = n, above = above, at_or_below = n - above, stringsAsFactors = FALSE)
  .rank_result("ad_median_test", "median", test, tab, pairs, groups, groups$median,
               paste0(if (exact) "Fisher exact test" else "median test", " by pairs (",
                      if (pair_median == "pair") "median of each pair" else "grand median", ")"),
               p_adjust, alpha, inp, list(grand_median = M, correct = correct, exact = exact))
}

#' @export
print.ad_median_test <- function(x, digits = 4, ...) {
  print.ad_nonparametric(x, digits = digits)
  cat(sprintf("\nGrand median = %s; counts in 'table' are above / at-or-below it.\n", format(x$grand_median, digits = digits)))
  if (any(x$table$n < 5)) cat("\nSome levels have fewer than 5 observations: chi-square p-values are rough; consider exact = TRUE.\n")
  invisible(x)
}
