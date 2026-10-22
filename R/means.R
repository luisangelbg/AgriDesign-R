# Mean separation: eleven pairwise procedures with compact letters, Dunnett
# against a control, simple effects in split- and strip-plots, orthogonal
# polynomial trends and user-defined contrasts.

.methods <- list(
  tukey = list(name = "Tukey HSD (Tukey-Kramer)", family = "q"),
  lsd = list(name = "Fisher's LSD", family = "t"),
  bonferroni = list(name = "Bonferroni", family = "t"),
  sidak = list(name = "Sidak", family = "t"),
  holm = list(name = "Holm (step-down Bonferroni)", family = "t"),
  duncan = list(name = "Duncan's multiple range test", family = "range"),
  snk = list(name = "Student-Newman-Keuls", family = "range"),
  regwq = list(name = "REGWQ (Ryan-Einot-Gabriel-Welsch)", family = "range"),
  scheffe = list(name = "Scheffe", family = "F"),
  dunnett = list(name = "Dunnett (vs control)", family = "dunnett"),
  games_howell = list(name = "Games-Howell", family = "gh"))

# ---------- Dunnett's multivariate t ----------
# Comparisons i = 1..k with the control have correlations lambda_i lambda_j,
# lambda_i = sqrt(n_i / (n_i + n_0)): one standard normal factor z makes them
# independent, T_i = (lambda_i z + sqrt(1 - lambda_i^2) e_i) / s, and
# s = sqrt(chi^2_df / df). Two-sided p = P(max |T_i| > t); the complement
# 1 - prod(...) is integrated directly (expm1 of a sum of logs) so that small
# p-values keep their relative precision. Both integrals use composite
# Gauss-Legendre rules (8 panels of 20 nodes), vectorised over z and s.

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
.gauss_legendre <- function(n) {
  b <- seq_len(n - 1) / sqrt(4 * seq_len(n - 1)^2 - 1)
  J <- matrix(0, n, n); J[cbind(1:(n - 1), 2:n)] <- b; J[cbind(2:n, 1:(n - 1))] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
}
.gl20 <- .gauss_legendre(20)

.composite_nodes <- function(lo, hi, panels = 8) {
  h <- (hi - lo) / panels; mid <- lo + h * (seq_len(panels) - 0.5)
  list(x = as.vector(outer(.gl20$x * h / 2, mid, `+`)), w = rep(.gl20$w * h / 2, panels))
}

.dunnett_p <- function(t, lambda, df) {
  t <- abs(t); s1 <- sqrt(1 - lambda^2)
  zn <- .composite_nodes(-10, 10)
  if (is.finite(df)) {
    sn <- .composite_nodes(max(0, 1 - 10 / sqrt(df)), 1 + 10 / sqrt(df))
    ws <- sn$w * exp(log(2) + (df / 2) * log(df / 2) - lgamma(df / 2) + (df - 1) * log(sn$x) - df * sn$x^2 / 2)
  } else {
    sn <- list(x = 1); ws <- 1
  }
  L <- matrix(0, length(zn$x), length(sn$x))
  for (i in seq_along(lambda)) {
    a <- outer(lambda[i] * zn$x, t * sn$x, `+`) / s1[i]
    b <- outer(lambda[i] * zn$x, -t * sn$x, `+`) / s1[i]
    L <- L + log(pmax(stats::pnorm(a) - stats::pnorm(b), 1e-300))
  }
  p <- sum((zn$w * stats::dnorm(zn$x)) * (-expm1(L)) %*% ws)
  min(1, max(0, p))
}

.dunnett_crit <- function(alpha, lambda, df) {
  stats::uniroot(function(t) .dunnett_p(t, lambda, df) - alpha, c(0.5, 10), tol = 1e-8)$root
}

# ---------- pairwise engine ----------
# g: data frame with label, mean, n, sd; se_diff(i, j); mse and dfe of the error used.
.compare <- function(g, se_diff, mse, dfe, alpha, method, control = 1) {
  k <- nrow(g)
  ij <- .pair_index(k)
  if (method == "dunnett") {
    others <- setdiff(seq_len(k), control)
    ij <- cbind(others, rep(control, length(others)))
  }
  m <- nrow(ij)
  diff <- g$mean[ij[, 1]] - g$mean[ij[, 2]]
  se <- vapply(seq_len(m), function(q) se_diff(ij[q, 1], ij[q, 2]), 0)
  tval <- diff / se; qval <- abs(diff) / se * sqrt(2)
  bal <- length(unique(g$n)) == 1
  se_bal <- if (bal) se_diff(1, 2) else NA_real_
  out <- data.frame(i = ij[, 1], j = ij[, 2], level_i = g$label[ij[, 1]], level_j = g$label[ij[, 2]],
                    diff = diff, se = se, stringsAsFactors = FALSE)
  msd <- NULL; note <- ""
  if (method %in% c("lsd", "bonferroni", "sidak", "holm")) {
    p <- 2 * stats::pt(abs(tval), dfe, lower.tail = FALSE)
    padj <- switch(method, lsd = p, bonferroni = pmin(1, p * m), sidak = 1 - (1 - p)^m,
                   holm = stats::p.adjust(p, "holm"))
    a <- switch(method, lsd = alpha, bonferroni = alpha / m, sidak = 1 - (1 - alpha)^(1 / m), holm = NA)
    out$statistic <- tval; out$p <- p; out$p_adj <- padj; out$different <- padj < alpha
    out$critical <- if (is.na(a)) NA_real_ else stats::qt(1 - a / 2, dfe) * se
    if (!is.na(a) && bal) msd <- c(value = stats::qt(1 - a / 2, dfe) * se_bal)
    if (method == "lsd") note <- "Per-comparison alpha: interpret as a protected LSD only if the F-test was significant."
  } else if (method == "tukey") {
    qc <- stats::qtukey(1 - alpha, k, dfe)
    p <- stats::ptukey(qval, k, dfe, lower.tail = FALSE)
    out$statistic <- qval; out$p <- p; out$p_adj <- p; out$different <- p < alpha
    out$critical <- qc * se / sqrt(2)
    if (bal) msd <- c(value = qc * se_bal / sqrt(2))
    if (!bal) note <- "Unequal replication: Tukey-Kramer form (own SE for each pair)."
  } else if (method == "scheffe") {
    Fs <- tval^2 / (k - 1)
    p <- stats::pf(Fs, k - 1, dfe, lower.tail = FALSE)
    cr <- sqrt((k - 1) * stats::qf(1 - alpha, k - 1, dfe))
    out$statistic <- Fs; out$p <- p; out$p_adj <- p; out$different <- p < alpha; out$critical <- cr * se
    if (bal) msd <- c(value = cr * se_bal)
  } else if (method %in% c("snk", "duncan", "regwq")) {
    ord <- order(-g$mean); pos <- integer(k); pos[ord] <- seq_len(k)
    alpha_p <- function(p) switch(method, snk = alpha, duncan = 1 - (1 - alpha)^(p - 1),
                                  regwq = if (p >= k - 1) alpha else 1 - (1 - alpha)^(p / k))
    homog <- matrix(FALSE, 0, 2)
    declared <- function(r1, r2) any(homog[, 1] <= r1 & r2 <= homog[, 2])
    qpair <- function(a, b) {
      w <- which((ij[, 1] == a & ij[, 2] == b) | (ij[, 1] == b & ij[, 2] == a))
      qval[w]
    }
    for (span in k:2) for (r1 in seq_len(k - span + 1)) {
      r2 <- r1 + span - 1
      if (declared(r1, r2)) next
      if (qpair(ord[r1], ord[r2]) < stats::qtukey(1 - alpha_p(span), span, dfe)) homog <- rbind(homog, c(r1, r2))
    }
    span <- abs(pos[ij[, 1]] - pos[ij[, 2]]) + 1
    ap <- vapply(span, alpha_p, 0)
    p <- stats::ptukey(qval, span, dfe, lower.tail = FALSE)
    out$statistic <- qval; out$span <- span; out$p <- p
    out$p_adj <- if (method == "duncan") 1 - (1 - p)^(1 / (span - 1)) else p
    out$different <- !mapply(function(a, b) declared(min(a, b), max(a, b)), pos[ij[, 1]], pos[ij[, 2]])
    out$critical <- stats::qtukey(1 - ap, span, dfe) * se / sqrt(2)
    if (bal) msd <- stats::setNames(vapply(2:k, function(p) stats::qtukey(1 - alpha_p(p), p, dfe) * se_bal / sqrt(2), 0),
                                    paste0("p=", 2:k))
    note <- switch(method,
      duncan = "Duncan uses the protection level (1 - alpha)^(p - 1) for a range of p means, which makes it liberal.",
      snk = "Step-down test: a difference inside a range already declared homogeneous is not tested.",
      regwq = "REGWQ uses alpha_p = 1 - (1 - alpha)^(p/k) for p < k - 1 and alpha for the two widest ranges.")
  } else if (method == "dunnett") {
    lambda <- sqrt(g$n[ij[, 1]] / (g$n[ij[, 1]] + g$n[control]))
    cr <- .dunnett_crit(alpha, lambda, dfe)
    p <- vapply(tval, .dunnett_p, 0, lambda = lambda, df = dfe)
    out$statistic <- tval; out$p <- p; out$p_adj <- p; out$different <- p < alpha; out$critical <- cr * se
    if (bal) msd <- c(value = cr * se_bal)
    note <- paste0("Two-sided comparisons with the control (", g$label[control], ").")
  } else if (method == "games_howell") {
    vi <- g$sd[ij[, 1]]^2 / g$n[ij[, 1]]; vj <- g$sd[ij[, 2]]^2 / g$n[ij[, 2]]
    out$se <- sqrt(vi + vj)
    df <- (vi + vj)^2 / (vi^2 / (g$n[ij[, 1]] - 1) + vj^2 / (g$n[ij[, 2]] - 1))
    qv <- abs(diff) / out$se * sqrt(2)
    p <- stats::ptukey(qv, k, df, lower.tail = FALSE)
    out$statistic <- qv; out$df <- df; out$p <- p; out$p_adj <- p; out$different <- p < alpha
    out$critical <- stats::qtukey(1 - alpha, k, df) * out$se / sqrt(2)
    note <- "Each pair uses its own SE and Welch degrees of freedom; no pooled error."
  }
  groups <- if (method == "dunnett") {
    ifelse(seq_len(k) == control, "control", ifelse(seq_len(k) %in% out$i[out$different], "*", "ns"))
  } else .letters(g$mean, .diff_matrix(k, ij, ifelse(out$different, 0, 1), 0.5))
  list(pairs = out, groups = groups, msd = msd, note = note)
}

.check_fit <- function(fit) if (!inherits(fit, "ad_anova")) stop("'fit' must come from design_anova()", call. = FALSE)

.check_factors <- function(fit, f) {
  facs <- fit$roles$treatment
  bad <- setdiff(f, facs)
  if (length(bad)) stop("not a treatment factor: ", paste(bad, collapse = ", "), call. = FALSE)
}

# SE of the difference of two means: covariance form under the residual error,
# balanced form MS (1/n_i + 1/n_j) under a higher stratum.
.se_diff_fun <- function(fit, lsm, err_ms, residual) {
  if (residual) {
    mse <- fit$residual[["ms"]]
    function(i, j) { dx <- lsm$xbar[i, ] - lsm$xbar[j, ]; sqrt(max(0, drop(t(dx) %*% lsm$XtXi %*% dx)) * mse) }
  } else {
    n <- lsm$means$n
    function(i, j) sqrt(err_ms * (1 / n[i] + 1 / n[j]))
  }
}

# Pooled error of simple effects in split- and strip-plots (Satterthwaite df).
.pooled <- function(main, sub, b) {
  ms <- ((b - 1) * sub$ms + main$ms) / b
  df <- ms^2 / (((b - 1) * sub$ms / b)^2 / sub$df + (main$ms / b)^2 / main$df)
  list(ms = ms, df = df)
}

#' Compare treatment means
#'
#' Separates the least-squares means of a treatment factor (or of the
#' combinations of two factors) with one of eleven procedures and labels them
#' with letters: levels that share a letter do not differ. Each factor is
#' compared with its own error stratum (for instance, main plots of a
#' split-plot with Error a).
#'
#' Methods: `"tukey"` (Tukey-Kramer HSD), `"lsd"` (Fisher), `"bonferroni"`,
#' `"sidak"`, `"holm"`, `"duncan"`, `"snk"` (Student-Newman-Keuls),
#' `"regwq"`, `"scheffe"`, `"dunnett"` (each level against `control`; the
#' groups are then `"*"` for a difference with the control and `"ns"`) and
#' `"games_howell"` (unequal variances; own SE and Welch df per pair).
#'
#' With `by`, the comparison is made within each level of another factor
#' (simple effects). In split-plots the sub-plot factor within a main plot
#' uses Error b, and the main-plot factor within a sub-plot level uses the
#' pooled error `[(b - 1) Eb + Ea] / b` with Satterthwaite degrees of
#' freedom; strip-plots pool their errors the same way.
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param factor Treatment factor whose means are compared; two names compare
#'   all their combinations together (not available for split- or
#'   strip-plots, where `by` should be used).
#' @param method Procedure, see Details.
#' @param alpha Significance level.
#' @param control Level used as control by Dunnett's test (default: first level).
#' @param by Optional second treatment factor: compare `factor` within each of
#'   its levels.
#' @return An object of class `"ad_comparison"` with the means and their
#'   letters (`means`), the pairwise comparisons (`pairs`), the minimum
#'   significant difference (`msd`, for equal replication) and the error used.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' compare_means(fit, "Nitrogen_kg_ha")
#'
#' fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
#' compare_means(fit, "Variety", method = "dunnett", control = "Negro Jamapa")
#'
#' fit <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
#'                     block = "Block", design = "split_rcbd")
#' compare_means(fit, "Variety", by = "Irrigation")
#' @export
compare_means <- function(fit, factor, method = "tukey", alpha = 0.05, control = NULL, by = NULL) {
  .check_fit(fit)
  method <- match.arg(method, names(.methods))
  .check_factors(fit, c(factor, by))
  if (!is.null(by)) return(.simple_effects(fit, factor, by, method, alpha, control))
  if (length(factor) > 2) stop("compare at most two factors at once; use 'by' for simple effects", call. = FALSE)
  split <- fit$design %in% c("split_rcbd", "split_crd", "strip_rcbd", "splitsplit_rcbd")
  if (length(factor) == 2 && split)
    stop("in split- and strip-plots compare the combinations with 'by' (simple effects)", call. = FALSE)
  lsm <- .lsm(fit, factor)
  err <- if (lsm$error == "Residuals") list(ms = fit$residual[["ms"]], df = fit$residual[["df"]], label = fit$strata$Residuals$label)
         else fit$strata[[lsm$error]]
  res <- .run_compare(fit, lsm, err, lsm$error == "Residuals", method, alpha, control)
  res$factor <- factor; res$by <- NULL
  res
}

.run_compare <- function(fit, lsm, err, residual, method, alpha, control, se_diff = NULL) {
  g <- lsm$means
  ctl <- 1
  if (method == "dunnett" && !is.null(control)) {
    ctl <- match(control, g$label)
    if (is.na(ctl)) stop("control level '", control, "' not found", call. = FALSE)
  }
  if (method == "games_howell" && any(g$n < 2)) stop("Games-Howell needs at least two observations per level", call. = FALSE)
  if (is.null(se_diff)) se_diff <- .se_diff_fun(fit, lsm, err$ms, residual)
  cmp <- .compare(g, se_diff, err$ms, err$df, alpha, method, ctl)
  g$group <- cmp$groups
  if (!is.null(fit$transform)) g$back_transformed <- fit$transform$inv(g$mean)
  structure(list(means = g, pairs = cmp$pairs, msd = cmp$msd, note = cmp$note,
                 method = method, method_name = .methods[[method]]$name, alpha = alpha,
                 error = list(ms = err$ms, df = err$df, label = err$label),
                 response = fit$response, transform = fit$transform),
            class = "ad_comparison")
}

.simple_effects <- function(fit, factor, by, method, alpha, control) {
  if (length(factor) != 1 || length(by) != 1 || factor == by)
    stop("'factor' and 'by' must be two different treatment factors", call. = FALSE)
  lsm <- .lsm(fit, c(factor, by))
  cells <- lsm$means
  tr <- fit$roles$treatment
  r <- fit$n / nrow(cells)
  spec <- NULL
  if (fit$design %in% c("split_rcbd", "split_crd") && all(c(factor, by) %in% tr[1:2])) {
    ea <- fit$strata[[fit$table$term[fit$table$is_error][1]]]
    eb <- fit$strata$Residuals
    if (factor == tr[2]) {
      spec <- list(ms = eb$ms, df = eb$df, label = "Error b (sub-plot)")
    } else {
      pe <- .pooled(ea, eb, nlevels(fit$data[[tr[2]]]))
      spec <- list(ms = pe$ms, df = pe$df,
                   label = sprintf("pooled error [(b-1)Eb + Ea]/b with Satterthwaite df"))
    }
  } else if (fit$design == "strip_rcbd") {
    errs <- fit$table$term[fit$table$is_error]
    ea <- fit$strata[[errs[1]]]; eb <- fit$strata[[errs[2]]]; ec <- fit$strata$Residuals
    if (factor == tr[2]) {
      pe <- .pooled(eb, ec, nlevels(fit$data[[tr[1]]]))
      spec <- list(ms = pe$ms, df = pe$df, label = sprintf("pooled error [(a-1)Ec + Eb]/a with Satterthwaite df"))
    } else {
      pe <- .pooled(ea, ec, nlevels(fit$data[[tr[2]]]))
      spec <- list(ms = pe$ms, df = pe$df, label = sprintf("pooled error [(b-1)Ec + Ea]/b with Satterthwaite df"))
    }
  }
  e_name <- .error_of(fit, c(factor, by))
  base_err <- if (e_name == "Residuals") list(ms = fit$residual[["ms"]], df = fit$residual[["df"]], label = "Residual error")
              else fit$strata[[e_name]]
  slices <- lapply(levels(fit$data[[by]]), function(bl) {
    sel <- which(cells[[by]] == bl)
    sub <- list(means = cells[sel, , drop = FALSE], xbar = lsm$xbar[sel, , drop = FALSE], XtXi = lsm$XtXi)
    sub$means$label <- as.character(sub$means[[factor]])
    err <- if (is.null(spec)) base_err else spec
    sd_fun <- if (is.null(spec)) .se_diff_fun(fit, sub, err$ms, e_name == "Residuals")
              else function(i, j) sqrt(2 * spec$ms / r)
    out <- .run_compare(fit, sub, err, e_name == "Residuals", method, alpha, control, se_diff = sd_fun)
    out$means[[by]] <- bl; out$pairs[[by]] <- rep(bl, nrow(out$pairs))
    out
  })
  first <- slices[[1]]
  structure(list(means = do.call(rbind, lapply(slices, `[[`, "means")),
                 pairs = do.call(rbind, lapply(slices, `[[`, "pairs")),
                 msd = first$msd, note = first$note, method = method, method_name = first$method_name,
                 alpha = alpha, error = first$error, response = fit$response, transform = fit$transform,
                 factor = factor, by = by),
            class = "ad_comparison")
}

#' @export
print.ad_comparison <- function(x, digits = 4, ...) {
  cat("\n", x$method_name, ", alpha = ", x$alpha, "\n", sep = "")
  cat("Response: ", x$response, "   compared: ", paste(x$factor, collapse = " x "),
      if (!is.null(x$by)) paste0(" within each level of ", x$by), "\n", sep = "")
  cat(sprintf("Error: %s (MS = %s, df = %s)\n", x$error$label, format(x$error$ms, digits = digits),
              format(x$error$df, digits = digits)))
  if (nzchar(x$note)) cat(x$note, "\n")
  if (!is.null(x$msd)) {
    if (length(x$msd) == 1) cat(sprintf("Minimum significant difference = %s\n", format(x$msd, digits = digits)))
    else cat("Least significant ranges: ", paste(names(x$msd), format(x$msd, digits = digits), sep = ": ", collapse = "; "), "\n")
  }
  m <- x$means
  cols <- c(if (!is.null(x$by)) x$by, "label", "n", "mean", "se", "raw_mean", "sd",
            if ("back_transformed" %in% names(m)) "back_transformed", "group")
  out <- m[cols]
  names(out)[names(out) == "label"] <- paste(x$factor, collapse = " x ")
  ord <- if (is.null(x$by)) order(-m$mean) else order(match(m[[x$by]], unique(m[[x$by]])), -m$mean)
  cat("\n"); print(out[ord, ], row.names = FALSE, digits = digits)
  if (x$method != "dunnett") cat("\nMeans sharing a letter do not differ", if (!is.null(x$by)) " (letters valid within each level)", ".\n", sep = "")
  invisible(x)
}

# ---------- contrasts and trends ----------

.contrast_one <- function(lsm, coef, err, residual) {
  est <- sum(coef * lsm$means$mean)
  v <- if (residual) { cx <- drop(coef %*% lsm$xbar); drop(t(cx) %*% lsm$XtXi %*% cx) }
       else sum(coef^2 / lsm$means$n)
  se <- sqrt(err$ms * v)
  tt <- est / se; ss <- est^2 / v
  c(estimate = est, se = se, t = tt, p = 2 * stats::pt(abs(tt), err$df, lower.tail = FALSE),
    ss = ss, F = ss / err$ms, df_error = err$df)
}

.factor_error <- function(fit, lsm) {
  if (lsm$error == "Residuals") list(ms = fit$residual[["ms"]], df = fit$residual[["df"]], label = "Residual error")
  else fit$strata[[lsm$error]]
}

#' Test contrasts among the levels of a factor
#'
#' A contrast compares groups of treatments with coefficients that sum to
#' zero, one per level in the order of `levels(fit$data[[factor]])`; for
#' instance a control against the other three treatments is `c(3, -1, -1, -1)`.
#'
#' @inheritParams compare_means
#' @param coef A vector of coefficients, or a matrix with one contrast per
#'   row (row names are used as labels).
#' @return A data frame with the estimate, its SE, t, p, the sum of squares
#'   of the contrast and F (1 and `df_error` degrees of freedom).
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' test_contrasts(fit, "Nitrogen_kg_ha", rbind("0 vs fertilised" = c(3, -1, -1, -1),
#'                                       "120 vs 180" = c(0, 0, 1, -1)))
#' @export
test_contrasts <- function(fit, factor, coef) {
  .check_fit(fit); .check_factors(fit, factor)
  if (length(factor) != 1) stop("contrasts are defined on one factor", call. = FALSE)
  lsm <- .lsm(fit, factor)
  k <- nrow(lsm$means)
  C <- if (is.matrix(coef)) coef else matrix(coef, nrow = 1)
  if (ncol(C) != k) stop("give ", k, " coefficients, one per level: ", paste(lsm$means$label, collapse = ", "), call. = FALSE)
  if (any(abs(rowSums(C)) > 1e-9)) stop("the coefficients of a contrast must sum to zero", call. = FALSE)
  err <- .factor_error(fit, lsm)
  res <- t(apply(C, 1, .contrast_one, lsm = lsm, err = err, residual = lsm$error == "Residuals"))
  data.frame(contrast = if (!is.null(rownames(C))) rownames(C) else apply(C, 1, paste, collapse = ", "),
             res, row.names = NULL, check.names = FALSE)
}

# Orthogonal polynomial coefficients for (possibly unequally spaced) levels
# with weights n, scaled to small integers when possible.
.ortho_poly <- function(x, n, max_deg) {
  w <- n; P <- list(rep(1, length(x)))
  dot <- function(a, b) sum(w * a * b)
  for (d in seq_len(max_deg)) {
    v <- x^d
    for (p in P) v <- v - dot(v, p) / dot(p, p) * p
    s <- v / min(abs(v[abs(v) > 1e-9]))
    for (mult in c(1, 2, 3, 4, 5, 6, 8, 10, 12)) {
      tt <- s * mult
      if (all(abs(tt - round(tt)) < 1e-6)) { s <- round(tt); break }
    }
    P[[d + 1]] <- s
  }
  P[-1]
}

#' Trend analysis of a quantitative factor
#'
#' Splits the sum of squares of a factor with numeric levels (doses, rates,
#' dates) into linear, quadratic, cubic... components with orthogonal
#' polynomial contrasts (also with unequal spacing) and fits the polynomial
#' of the highest significant degree to the observations. A fitted quadratic
#' reports where the maximum (or minimum) response lies.
#'
#' @inheritParams compare_means
#' @param degree Highest degree to test (default: up to 4, at most the
#'   number of levels minus one).
#' @return An object of class `"ad_trend"` with the table of components
#'   (`components`), the fitted coefficients (`coef`, intercept first), its
#'   R squared and the optimum of a quadratic (`optimum`).
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' trend(fit, "Nitrogen_kg_ha")
#' @export
trend <- function(fit, factor, degree = NULL, alpha = 0.05) {
  .check_fit(fit); .check_factors(fit, factor)
  lv <- levels(fit$data[[factor]])
  x <- suppressWarnings(as.numeric(lv))
  if (length(lv) < 3 || any(is.na(x))) stop("trend analysis needs at least three numeric levels", call. = FALSE)
  max_deg <- min(if (is.null(degree)) 4 else degree, length(lv) - 1)
  lsm <- .lsm(fit, factor)
  err <- .factor_error(fit, lsm)
  P <- .ortho_poly(x, lsm$means$n, max_deg)
  nm <- c("Linear", "Quadratic", "Cubic", "Quartic", "Quintic")
  comp <- do.call(rbind, lapply(seq_along(P), function(d) {
    r <- .contrast_one(lsm, P[[d]], err, lsm$error == "Residuals")
    data.frame(component = if (d <= 5) nm[d] else paste("Degree", d), coefficients = paste(P[[d]], collapse = ", "),
               ss = r[["ss"]], F = r[["F"]], df_error = r[["df_error"]],
               p = stats::pf(r[["F"]], 1, r[["df_error"]], lower.tail = FALSE), stringsAsFactors = FALSE)
  }))
  deg <- if (any(comp$p < alpha)) max(which(comp$p < alpha)) else 0
  fit_deg <- max(1, deg)
  xs <- as.numeric(as.character(fit$data[[factor]]))
  pf <- stats::lm.fit(outer(xs, 0:fit_deg, `^`), fit$y)
  b <- unname(pf$coefficients)
  r2 <- 1 - sum(pf$residuals^2) / sum((fit$y - mean(fit$y))^2)
  opt <- NULL
  if (fit_deg == 2 && b[3] != 0) {
    xo <- -b[2] / (2 * b[3])
    opt <- c(x = xo, y = b[1] + b[2] * xo + b[3] * xo^2)
    attr(opt, "kind") <- if (b[3] < 0) "maximum" else "minimum"
  }
  structure(list(components = comp, significant_degree = deg, coef = b, r_squared = r2, optimum = opt,
                 range = range(x), factor = factor, response = fit$response, error = err, alpha = alpha),
            class = "ad_trend")
}

#' @export
print.ad_trend <- function(x, digits = 4, ...) {
  cat("\nTrend analysis of ", x$response, " on ", x$factor, " (error: ", x$error$label, ")\n\n", sep = "")
  print(x$components, row.names = FALSE, digits = digits)
  b <- x$coef
  terms <- c(format(b[1], digits = digits), vapply(seq_along(b)[-1], function(i)
    paste0(if (b[i] < 0) "- " else "+ ", format(abs(b[i]), digits = digits), " x", if (i > 2) paste0("^", i - 1) else ""), ""))
  cat("\n", if (x$significant_degree) paste0("Highest significant component: ", tolower(x$components$component[x$significant_degree]), ".")
      else "No polynomial component is significant.", "\n", sep = "")
  cat("Fitted response: y = ", paste(terms, collapse = " "), sprintf("   R-squared = %.3f\n", x$r_squared), sep = "")
  if (!is.null(x$optimum)) {
    o <- x$optimum
    kind <- attr(o, "kind")
    cat(sprintf("%s response at x = %s (y = %s)%s\n", paste0(toupper(substr(kind, 1, 1)), substring(kind, 2)),
                format(o[["x"]], digits = digits),
                format(o[["y"]], digits = digits),
                if (o[["x"]] < x$range[1] || o[["x"]] > x$range[2]) "; outside the levels tested, do not extrapolate"
                else if (o[["x"]] > x$range[2] - 0.05 * diff(x$range)) "; at the edge of the levels tested" else ""))
  }
  invisible(x)
}
