# Assumptions of the analysis of variance: normality of the residuals,
# homogeneity of variances among treatments, additivity of blocks and
# treatments, independence and extreme observations; transformations and
# the Box-Cox profile.

# ---------- tests that base R does not provide ----------

.levene <- function(x, g, center = c("median", "mean")) {
  center <- match.arg(center)
  cfun <- if (center == "median") stats::median else mean
  z <- abs(x - stats::ave(x, g, FUN = cfun))
  r <- .oneway(z, g)
  c(F = r[["F"]], df1 = r[["df1"]], df2 = r[["df2"]], p = r[["p"]])
}

.oneway <- function(x, g) {
  g <- factor(g)
  m <- tapply(x, g, mean); n <- tabulate(g)
  ssb <- sum(n * (m - mean(x))^2); ssw <- sum((x - m[g])^2)
  df1 <- nlevels(g) - 1; df2 <- length(x) - nlevels(g)
  F <- (ssb / df1) / (ssw / df2)
  c(F = F, df1 = df1, df2 = df2, p = stats::pf(F, df1, df2, lower.tail = FALSE))
}

# Anderson-Darling with estimated mean and variance; Stephens' p for the
# modified statistic. Logs of the normal tails are taken without losing them.
.anderson_darling <- function(x) {
  x <- sort(x); n <- length(x)
  if (n < 8) return(c(A = NA_real_, p = NA_real_))
  z <- (x - mean(x)) / stats::sd(x)
  lp <- stats::pnorm(z, log.p = TRUE)
  lq <- stats::pnorm(rev(z), lower.tail = FALSE, log.p = TRUE)
  A <- -n - sum((2 * seq_len(n) - 1) * (lp + lq)) / n
  AA <- A * (1 + 0.75 / n + 2.25 / n^2)
  p <- if (AA < 0.2) 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
       else if (AA < 0.34) 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
       else if (AA < 0.6) exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
       else if (AA < 10) exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
       else 3.7e-24
  c(A = A, p = min(1, max(0, p)))
}

.jarque_bera <- function(x) {
  n <- length(x); d <- x - mean(x)
  m2 <- mean(d^2); m3 <- mean(d^3); m4 <- mean(d^4)
  sk <- m3 / m2^1.5; ku <- m4 / m2^2 - 3
  JB <- n / 6 * (sk^2 + ku^2 / 4)
  c(JB = JB, p = stats::pchisq(JB, 2, lower.tail = FALSE), skewness = sk, kurtosis = ku)
}

# Sample skewness with the small-sample factor n / ((n - 1)(n - 2)).
.skewness <- function(x) {
  n <- length(x); s <- stats::sd(x)
  if (n < 3) return(NA_real_)
  if (!s) return(0)
  n / ((n - 1) * (n - 2)) * sum(((x - mean(x)) / s)^3)
}

# Tukey's one-degree-of-freedom test: adds (fitted - mean)^2 to the model.
.tukey_additivity <- function(X, y) {
  f1 <- .ls(X, y)
  if (f1$df < 2) return(NULL)
  fv <- y - f1$fit$residuals
  f2 <- .ls(cbind(X, (fv - mean(fv))^2), y)
  df <- f1$df - f2$df
  if (df < 1) return(NULL)
  ss <- f1$sse - f2$sse
  F <- (ss / df) / (f2$sse / f2$df)
  c(F = F, df1 = df, df2 = f2$df, p = stats::pf(F, df, f2$df, lower.tail = FALSE), ss = ss)
}

.durbin_watson <- function(e) sum(diff(e)^2) / sum(e^2)

.hat_values <- function(f) {
  Q <- qr.Q(f$qr)[, seq_len(f$rank), drop = FALSE]
  rowSums(Q^2)
}

# Traffic light of a p-value: ok, warn (between alpha/5 and alpha) or bad.
.level <- function(p, alpha) {
  if (is.null(p) || !is.finite(p)) "info" else if (p >= alpha) "ok" else if (p >= alpha / 5) "warn" else "bad"
}

.treatment_groups <- function(fit) {
  tr <- fit$roles$treatment
  if (length(tr) == 1) fit$data[[tr]] else interaction(fit$data[tr], sep = " x ", drop = TRUE)
}

.has_blocking <- function(fit) length(fit$roles$block) > 0 || (!is.null(fit$roles$row) && !is.null(fit$roles$col))

# Every test on one response vector, with the model matrix of the design.
.diagnose <- function(X, y, groups, blocking) {
  f <- .ls(X, y)
  e <- f$fit$residuals
  mse <- f$sse / f$df
  h <- .hat_values(f$fit)
  stud <- e / (sqrt(mse) * sqrt(pmax(1e-12, 1 - h)))
  can_var <- nlevels(groups) > 1 && all(table(groups) > 1)
  list(
    fit = f, resid = e, studentized = stud, mse = mse,
    shapiro = if (length(e) >= 3 && length(e) <= 5000) stats::shapiro.test(e) else NULL,
    anderson = .anderson_darling(e),
    jarque = .jarque_bera(e),
    levene = if (can_var) .levene(e, groups, "median") else NULL,
    levene_mean = if (can_var) .levene(e, groups, "mean") else NULL,
    bartlett = if (can_var) stats::bartlett.test(y, groups) else NULL,
    fligner = if (can_var) stats::fligner.test(y, groups) else NULL,
    fmax = if (can_var) { v <- tapply(y, groups, stats::var); max(v) / min(v) } else NULL,
    additivity = if (blocking) .tukey_additivity(X, y) else NULL,
    dw = .durbin_watson(e),
    cv = sqrt(mse) / abs(mean(y)) * 100)
}

#' Check the assumptions of the analysis of variance
#'
#' Runs the usual checks on the residuals of a fitted design: normality
#' (Shapiro-Wilk, Anderson-Darling, Jarque-Bera), homogeneity of variances
#' among treatments (Levene with median and with mean, Bartlett,
#' Fligner-Killeen, Hartley's Fmax), additivity of blocks and treatments
#' (Tukey's one-degree-of-freedom test, only for blocked designs),
#' independence in data order (Durbin-Watson) and extreme studentized
#' residuals. The checks are summarised in a grade from A (assumptions
#' satisfied) to D (several assumptions fail).
#'
#' Each check is rated "ok" when its p-value is at least `alpha`, "warn"
#' when it lies between `alpha / 5` and `alpha` and "bad" below `alpha / 5`.
#' The grade is D with two or more "bad" checks, C with one "bad" or two
#' "warn", B with a single "warn" and A otherwise.
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param alpha Significance level of the checks.
#' @return An object of class `"ad_assumptions"` with the table of tests
#'   (`tests`), the residuals and studentized residuals (`residuals`), the
#'   observations with |studentized residual| > 2.5 (`suspicious`) and the
#'   grade.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' check_assumptions(fit)
#' @export
check_assumptions <- function(fit, alpha = 0.05) {
  if (!inherits(fit, "ad_anova")) stop("'fit' must come from design_anova()", call. = FALSE)
  groups <- .treatment_groups(fit)
  dg <- .diagnose(fit$X, fit$y, groups, .has_blocking(fit))
  row <- function(check, test, stat_name, stat, df, p, level) {
    data.frame(check = check, test = test, statistic_name = stat_name, statistic = stat,
               df = df, p = p, level = level, stringsAsFactors = FALSE)
  }
  tests <- rbind(
    row("normality", "Shapiro-Wilk", "W",
        if (is.null(dg$shapiro)) NA else unname(dg$shapiro$statistic), "",
        if (is.null(dg$shapiro)) NA else dg$shapiro$p.value,
        .level(if (is.null(dg$shapiro)) NA else dg$shapiro$p.value, alpha)),
    row("normality", "Anderson-Darling", "A2", dg$anderson[["A"]], "", dg$anderson[["p"]], "info"),
    row("normality", "Jarque-Bera", "JB", dg$jarque[["JB"]], "2", dg$jarque[["p"]], "info"))
  if (!is.null(dg$levene)) {
    tests <- rbind(tests,
      row("variances", "Levene (median)", "F", dg$levene[["F"]],
          sprintf("%d, %d", dg$levene[["df1"]], dg$levene[["df2"]]), dg$levene[["p"]], .level(dg$levene[["p"]], alpha)),
      row("variances", "Levene (mean)", "F", dg$levene_mean[["F"]],
          sprintf("%d, %d", dg$levene_mean[["df1"]], dg$levene_mean[["df2"]]), dg$levene_mean[["p"]], "info"),
      row("variances", "Bartlett", "K2", unname(dg$bartlett$statistic),
          as.character(dg$bartlett$parameter), dg$bartlett$p.value, "info"),
      row("variances", "Fligner-Killeen", "X2", unname(dg$fligner$statistic),
          as.character(dg$fligner$parameter), dg$fligner$p.value, "info"),
      row("variances", "Hartley Fmax", "Fmax", dg$fmax, "", NA, "info"))
  }
  if (!is.null(dg$additivity)) {
    a <- dg$additivity
    tests <- rbind(tests, row("additivity", "Tukey one-degree-of-freedom", "F", a[["F"]],
                              sprintf("%d, %d", a[["df1"]], a[["df2"]]), a[["p"]], .level(a[["p"]], alpha)))
  }
  big <- sum(abs(dg$studentized) > 3)
  tests <- rbind(tests,
    row("extremes", "|studentized residual| > 3", "count", big, "", NA, if (big) "warn" else "ok"),
    row("independence", "Durbin-Watson (data order)", "d", dg$dw, "", NA,
        if (dg$dw > 1.5 && dg$dw < 2.5) "ok" else "info"))
  rated <- tests$level[tests$test %in% c("Shapiro-Wilk", "Levene (median)", "Tukey one-degree-of-freedom",
                                         "|studentized residual| > 3", "Durbin-Watson (data order)")]
  n_bad <- sum(rated == "bad"); n_warn <- sum(rated == "warn")
  grade <- if (n_bad >= 2) "D" else if (n_bad == 1 || n_warn >= 2) "C" else if (n_warn == 1) "B" else "A"
  res <- data.frame(row = as.integer(rownames(fit$data)), group = as.character(groups),
                    observed = fit$y, fitted = fit$y - dg$resid, residual = dg$resid,
                    studentized = dg$studentized, stringsAsFactors = FALSE)
  structure(list(tests = tests, grade = grade, alpha = alpha, residuals = res,
                 suspicious = res[abs(res$studentized) > 2.5, , drop = FALSE],
                 skewness = dg$jarque[["skewness"]], kurtosis = dg$jarque[["kurtosis"]],
                 response = fit$response, has_variances = !is.null(dg$levene),
                 has_additivity = !is.null(dg$additivity)),
            class = "ad_assumptions")
}

.grade_text <- c(
  A = "Assumptions satisfied. Proceed with the analysis of variance on this scale.",
  B = "Acceptable: minor deviations. Mention the check in the methods and, if the deviation is in the variances, be careful with the mean comparisons.",
  C = "Problematic: clear deviation from at least one assumption. Compare the transformations (transformations()) or use the rank-based route (nonparametric()).",
  D = "Assumptions fail: several checks fail. Rank-based procedures (nonparametric()) are the safer choice, or a generalised linear model for counts and proportions.")

#' @export
print.ad_assumptions <- function(x, digits = 4, ...) {
  cat("\nAssumptions of the ANOVA for ", x$response, " (alpha = ", x$alpha, ")\n\n", sep = "")
  t <- x$tests
  mark <- c(ok = "ok", warn = "WARN", bad = "BAD", info = "")
  out <- data.frame(Check = t$check, Test = t$test,
                    Statistic = paste(t$statistic_name, "=", trimws(formatC(t$statistic, digits = digits, format = "g"))),
                    df = t$df, p = ifelse(is.na(t$p), "", format.pval(t$p, digits = digits)),
                    Rating = mark[t$level], stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = FALSE)
  if (!x$has_variances) cat("\nVariances not testable: each treatment has a single observation.\n")
  if (nrow(x$suspicious)) {
    cat("\nObservations with |studentized residual| > 2.5:\n")
    print(x$suspicious[c("row", "group", "observed", "fitted", "studentized")], row.names = FALSE, digits = digits)
  }
  cat("\nGrade ", x$grade, ". ", .grade_text[[x$grade]], "\n", sep = "")
  invisible(x)
}

# ---------- transformations ----------

.transform_list <- function(y) {
  mn <- min(y); mx <- max(y)
  L <- list(none = list(label = "None (original scale)", f = identity, inv = identity, when = ""))
  if (mn > 0) L$log <- list(label = "Natural logarithm ln(y)", f = log, inv = exp,
                            when = "multiplicative effects, CV constant across treatments, SD proportional to the mean")
  if (mn >= 0) L$log1 <- list(label = "ln(y + 1)", f = function(v) log(v + 1), inv = function(v) exp(v) - 1,
                              when = "counts with zeros")
  if (mn > 0) L$log10 <- list(label = "log10(y)", f = log10, inv = function(v) 10^v, when = "")
  if (mn >= 0) L$sqrt <- list(label = "Square root sqrt(y)", f = sqrt, inv = function(v) v^2,
                              when = "Poisson counts, variance proportional to the mean")
  if (mn >= 0) L$sqrt05 <- list(label = "sqrt(y + 0.5)", f = function(v) sqrt(v + 0.5), inv = function(v) v^2 - 0.5,
                                when = "small counts with zeros")
  if (mn >= 0 && mx <= 100 && mx > 1)
    L$asin <- list(label = "Arcsine sqrt(y/100)", f = function(v) asin(sqrt(v / 100)), inv = function(v) 100 * sin(v)^2,
                   when = "percentages from counts, especially below 20 % or above 80 %")
  if (mn >= 0 && mx <= 1)
    L$asinp <- list(label = "Arcsine sqrt(y) (proportion)", f = function(v) asin(sqrt(v)), inv = function(v) sin(v)^2,
                    when = "proportions from counts")
  if (mn > 0 && mx < 1)
    L$logit <- list(label = "Logit ln(y / (1 - y))", f = function(v) log(v / (1 - v)), inv = function(v) 1 / (1 + exp(-v)),
                    when = "proportions away from 0 and 1")
  if (mn > 0) L$inv <- list(label = "Reciprocal 1/y", f = function(v) 1 / v, inv = function(v) 1 / v,
                            when = "rates and times; SD proportional to the square of the mean")
  L
}

#' Box-Cox profile likelihood of a design
#'
#' Profile log-likelihood of the Box-Cox power `lambda` for the model of the
#' design, its maximum, the 95 % interval and the nearest conventional power
#' (-1, -0.5, 0, 0.5, 1 or 2).
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param lambda Grid of powers.
#' @return A list with the profile (`profile`), the estimate (`lambda`), the
#'   95 % interval (`lower`, `upper`), the nearest conventional power
#'   (`rounded`) and whether the interval contains 1 (`contains_one`); `NULL`
#'   when the response has zeros or negative values.
#' @examples
#' fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
#'                     block = "Block")
#' box_cox(fit)[c("lambda", "lower", "upper", "rounded")]
#' @export
box_cox <- function(fit, lambda = seq(-2, 2, by = 0.1)) {
  if (!inherits(fit, "ad_anova")) stop("'fit' must come from design_anova()", call. = FALSE)
  y <- fit$y
  if (min(y) <= 0) return(NULL)
  n <- length(y); slog <- sum(log(y))
  ll <- vapply(lambda, function(l) {
    yt <- if (abs(l) < 1e-9) log(y) else (y^l - 1) / l
    -n / 2 * log(.ls(fit$X, yt)$sse / n) + (l - 1) * slog
  }, 0)
  best <- which.max(ll)
  cut <- ll[best] - stats::qchisq(0.95, 1) / 2
  inside <- lambda[ll >= cut]
  cand <- c(-1, -0.5, 0, 0.5, 1, 2)
  rounded <- cand[which.min(abs(cand - lambda[best]))]
  list(profile = data.frame(lambda = lambda, loglik = ll), lambda = lambda[best],
       lower = min(inside), upper = max(inside), rounded = rounded, cut = cut,
       contains_one = min(inside) <= 1 && max(inside) >= 1)
}

#' Compare transformations of the response
#'
#' Refits the design after each transformation that suits the range of the
#' response and repeats the assumption checks. A transformation is suggested
#' only when the original scale fails some check and the transformed scale
#' passes all of them: the conventional power nearest to the Box-Cox
#' estimate when it passes, otherwise the candidate with the least evidence
#' against the assumptions.
#'
#' @inheritParams check_assumptions
#' @return An object of class `"ad_transformations"`: the table of candidates
#'   (`table`), the identifier of the suggested transformation (`suggested`,
#'   `NULL` when none is needed or none works), a short `verdict` and the
#'   Box-Cox profile (`box_cox`).
#' @examples
#' fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
#'                     block = "Block")
#' transformations(fit)
#' @export
transformations <- function(fit, alpha = 0.05) {
  if (!inherits(fit, "ad_anova")) stop("'fit' must come from design_anova()", call. = FALSE)
  groups <- .treatment_groups(fit); blocking <- .has_blocking(fit)
  L <- .transform_list(fit$y)
  rows <- lapply(names(L), function(id) {
    yt <- L[[id]]$f(fit$y)
    if (any(!is.finite(yt))) return(NULL)
    dg <- .diagnose(fit$X, yt, groups, blocking)
    sw_p <- if (is.null(dg$shapiro)) NA_real_ else dg$shapiro$p.value
    lev_p <- if (is.null(dg$levene)) NA_real_ else dg$levene[["p"]]
    add_p <- if (is.null(dg$additivity)) NA_real_ else dg$additivity[["p"]]
    data.frame(id = id, transformation = L[[id]]$label,
               shapiro_W = if (is.null(dg$shapiro)) NA_real_ else unname(dg$shapiro$statistic),
               p_normality = sw_p, p_levene = lev_p, p_additivity = add_p,
               fmax = if (is.null(dg$fmax)) NA_real_ else dg$fmax,
               skewness = .skewness(dg$resid), cv = dg$cv,
               score = min(c(sw_p, lev_p, add_p, 1), na.rm = TRUE),
               typical_use = L[[id]]$when, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  bc <- box_cox(fit)
  bc_id <- if (is.null(bc)) NA_character_ else c(`-1` = "inv", `0` = "log", `0.5` = "sqrt")[as.character(bc$rounded)]
  none <- tab[tab$id == "none", ]
  bc_row <- if (!is.na(bc_id) && bc_id %in% tab$id) tab[tab$id == bc_id, ] else NULL
  best <- if (!is.null(bc_row) && bc_row$score >= alpha) bc_row else tab[which.max(tab$score), ]
  worth <- best$id != "none" && best$score >= alpha && none$score < alpha
  verdict <- if (worth) {
    paste0("Suggested: ", best$transformation, ". ",
           if (!is.null(bc_row) && best$id == bc_row$id)
             "It is the power closest to the Box-Cox estimate and it satisfies the checks"
           else "It is the transformation with the least evidence against the assumptions",
           sprintf(" (smallest p = %.4g vs %.4g on the original scale).", best$score, none$score),
           " Means are compared on the transformed scale; report them back-transformed.")
  } else if (none$score >= alpha) {
    "No transformation needed: the original scale already satisfies the assumptions."
  } else {
    sprintf("No transformation fixes all the problems; the best candidate is %s (smallest p = %.4g). Consider nonparametric(), or a generalised linear model for counts and proportions.",
            best$transformation, best$score)
  }
  structure(list(table = tab, suggested = if (worth) best$id else NULL, verdict = verdict,
                 box_cox = bc, alpha = alpha, response = fit$response),
            class = "ad_transformations")
}

#' @export
print.ad_transformations <- function(x, digits = 4, ...) {
  cat("\nTransformations of ", x$response, " (same model, alpha = ", x$alpha, ")\n\n", sep = "")
  t <- x$table
  star <- ifelse(!is.null(x$suggested) & t$id %in% x$suggested, "* ", "  ")
  fp <- function(p) ifelse(is.na(p), "-", format.pval(p, digits = 2, eps = 1e-4))
  out <- data.frame(` ` = star, Transformation = t$transformation, W = round(t$shapiro_W, 3),
                    `p norm.` = fp(t$p_normality), `p Levene` = fp(t$p_levene),
                    `p addit.` = fp(t$p_additivity), Fmax = round(t$fmax, 2),
                    Skew = round(t$skewness, 2), `CV %` = round(t$cv, 1),
                    check.names = FALSE, stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = FALSE)
  cat("(p norm.: Shapiro-Wilk; p Levene: Levene with median; p addit.: Tukey's additivity test)\n")
  if (!is.null(x$box_cox)) {
    b <- x$box_cox
    cat(sprintf("\nBox-Cox: lambda = %.2f, 95 %% interval %.1f to %.1f; %s\n", b$lambda, b$lower, b$upper,
                if (b$contains_one) "the interval includes 1 (no transformation required by this criterion)."
                else paste0("nearest conventional power ", b$rounded, ".")))
  } else cat("\nBox-Cox needs strictly positive values.\n")
  cat("\n", x$verdict, "\n", sep = "")
  invisible(x)
}
