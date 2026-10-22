# Analysis of variance for the designs of the catalogue.
# Factors use sum-to-zero (effect) coding, so dropping one term from the full
# model gives Type III sums of squares; aliased columns are dropped by the
# pivoting QR of lm.fit(). Each term is then tested against its own error stratum.

# Level order: factor levels if the column is a factor, sorted values if it is
# numeric, order of appearance otherwise.
.levels_of <- function(x) {
  if (is.factor(x)) return(levels(droplevels(x)))
  if (is.numeric(x)) return(as.character(sort(unique(x))))
  unique(as.character(x))
}

.check_columns <- function(data, cols) {
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop("column(s) not found in the data: ", paste(miss, collapse = ", "), call. = FALSE)
}

# Complete cases of the response, the factors and the covariates; role columns
# become factors and covariates are centred (least-squares means are then
# computed at the mean of each covariate).
.model_frame <- function(data, response, r) {
  fac <- unique(c(r$block, r$row, r$col, r$treatment))
  cols <- c(response, fac, r$covariate)
  .check_columns(data, cols)
  d <- as.data.frame(data)[cols]
  y <- d[[response]]
  if (!is.numeric(y)) stop("the response '", response, "' must be numeric", call. = FALSE)
  ok <- stats::complete.cases(d)
  d <- d[ok, , drop = FALSE]
  for (f in fac) d[[f]] <- factor(as.character(d[[f]]), levels = .levels_of(data[[f]][ok]))
  for (x in r$covariate) {
    if (!is.numeric(d[[x]])) stop("the covariate '", x, "' must be numeric", call. = FALSE)
    d[[x]] <- d[[x]] - mean(d[[x]])
  }
  list(data = d, dropped = sum(!ok), factors = fac)
}

.formula_of <- function(response, terms) {
  rhs <- vapply(terms, function(t) paste(c(t$factors, t$numeric), collapse = ":"), "")
  stats::as.formula(paste(sprintf("`%s`", response), "~",
                          paste(gsub("([^:]+)", "`\\1`", rhs), collapse = " + ")))
}

.contrasts_for <- function(fac) stats::setNames(rep(list("contr.sum"), length(fac)), fac)

# Least squares on the columns `keep` of X.
.ls <- function(X, y, keep = seq_len(ncol(X))) {
  f <- stats::lm.fit(X[, keep, drop = FALSE], y)
  sse <- sum(f$residuals^2)
  list(rank = f$rank, df = length(y) - f$rank, sse = sse, fit = f, keep = keep)
}

#' Analysis of variance of an experimental design
#'
#' Fits the linear model of a classical agricultural design and tests every
#' term against the right error stratum (for instance, main plots of a
#' split-plot against Error a).
#'
#' @param data A data frame in long format: one row per plot.
#' @param response Name of the numeric response column.
#' @param treatment Name(s) of the treatment factor(s). In split-plots,
#'   strip-plots and nested designs the order matters: the first factor goes
#'   on main plots (or is the outer factor).
#' @param block Name of the block column; two names (replicate, block within
#'   replicate) for incomplete block designs.
#' @param row,col Names of the row and column factors of a Latin square.
#' @param covariate Name(s) of numeric covariates (analysis of covariance).
#' @param design Identifier of the design (see [designs()]). When `NULL` the
#'   design is suggested from the roles given.
#' @param ss_type Type of sums of squares: `3` (default, adjusted for all
#'   other terms) or `1` (sequential).
#' @param random_blocks In a generalised RCBD, treat blocks as random and test
#'   treatments against the block x treatment interaction.
#' @param random_nested In a nested design, treat the inner factor as random
#'   and test the outer factor against it.
#' @return An object of class `"ad_anova"`: a list with the ANOVA table
#'   (`table`), the error strata (`strata`), the grand mean, the coefficient
#'   of variation, R squared, the fitted model and the data used.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha",
#'                     block = "Block", design = "rcbd")
#' fit
#'
#' design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
#'              block = "Block", design = "split_rcbd")
#' @export
design_anova <- function(data, response, treatment, block = NULL, row = NULL, col = NULL,
                         covariate = NULL, design = NULL, ss_type = 3,
                         random_blocks = TRUE, random_nested = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame", call. = FALSE)
  if (!ss_type %in% c(1, 3)) stop("'ss_type' must be 1 or 3", call. = FALSE)
  r <- list(treatment = treatment, block = block, row = row, col = col, covariate = covariate)
  mf <- .model_frame(data, response, r)
  d <- mf$data
  app <- .applicable(r, d)
  if (is.null(design)) design <- app$suggested
  if (!design %in% names(.catalog)) stop("unknown design '", design, "'; see designs()", call. = FALSE)
  chk <- app$table[app$table$id == design, ]
  if (!chk$ok) stop("design '", design, "' does not fit these roles: ", chk$why, call. = FALSE)
  for (f in mf$factors) if (nlevels(d[[f]]) < 2) stop("factor '", f, "' has fewer than two levels", call. = FALSE)

  dsg <- .catalog[[design]]
  terms <- dsg$build(r, list(random_blocks = random_blocks, random_nested = random_nested))
  fml <- .formula_of(response, terms)
  # keep.order: the columns of X must follow the catalogue terms (Type I and strata)
  tt <- stats::terms(fml, keep.order = TRUE)
  X <- stats::model.matrix(tt, d, contrasts.arg = .contrasts_for(mf$factors))
  assign <- attr(X, "assign")
  y <- d[[response]]
  full <- .ls(X, y)
  if (full$df < 1) stop("no degrees of freedom left for the residual error", call. = FALSE)

  ss <- df <- numeric(length(terms))
  if (ss_type == 1) {
    prev <- .ls(X, y, which(assign == 0))
    for (i in seq_along(terms)) {
      cur <- .ls(X, y, which(assign <= i))
      ss[i] <- prev$sse - cur$sse; df[i] <- prev$df - cur$df
      prev <- cur
    }
  } else {
    for (i in seq_along(terms)) {
      red <- .ls(X, y, which(assign != i))
      ss[i] <- red$sse - full$sse; df[i] <- red$df - full$df
    }
  }
  ms <- ifelse(df > 0, ss / df, NA_real_)
  mse <- full$sse / full$df
  tab <- data.frame(
    term = vapply(terms, `[[`, "", "name"),
    source = vapply(terms, `[[`, "", "label"),
    df = df, ss = ss, ms = ms,
    error = vapply(terms, `[[`, "", "error"),
    is_error = vapply(terms, `[[`, NA, "is_error"),
    stringsAsFactors = FALSE)
  err_ms <- function(e) if (e == "Residuals") c(mse, full$df) else unlist(tab[tab$term == e, c("ms", "df")])
  em <- t(vapply(tab$error, err_ms, c(0, 0)))
  tab$ms_error <- em[, 1]; tab$df_error <- em[, 2]
  tab$F <- ifelse(tab$df > 0 & tab$df_error > 0, tab$ms / tab$ms_error, NA_real_)
  tab$p <- stats::pf(tab$F, tab$df, tab$df_error, lower.tail = FALSE)
  # error strata are not tested themselves
  tab$F[tab$is_error] <- NA_real_; tab$p[tab$is_error] <- NA_real_

  strata <- list()
  for (i in which(tab$is_error)) strata[[tab$term[i]]] <- list(ms = tab$ms[i], df = tab$df[i], label = tab$source[i])
  strata$Residuals <- list(ms = mse, df = full$df,
                           label = if (length(strata)) "Residual error (last stratum)" else "Residual error")
  gm <- mean(y)
  sst <- sum((y - gm)^2)
  structure(list(
    table = tab, residual = c(df = full$df, ss = full$sse, ms = mse),
    total = c(df = length(y) - 1, ss = sst), strata = strata,
    design = design, design_name = dsg$name, response = response, roles = r,
    terms = terms, formula = fml, model_terms = tt, ss_type = ss_type,
    grand_mean = gm, cv = sqrt(mse) / gm * 100, r_squared = 1 - full$sse / sst,
    n = length(y), dropped = mf$dropped,
    X = X, y = y, qr_fit = full$fit, data = d),
    class = "ad_anova")
}

.stars <- function(p) {
  s <- stats::symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
              symbols = c("***", "**", "*", ".", " "))
  ifelse(is.na(p), "", as.character(s))
}

#' @export
print.ad_anova <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n", x$design_name, "\n", sep = "")
  cat("Response: ", x$response, "   (n = ", x$n,
      if (x$dropped) paste0(", ", x$dropped, " incomplete row(s) dropped"), ")\n", sep = "")
  cat("Type ", if (x$ss_type == 1) "I" else "III", " sums of squares\n\n", sep = "")
  t <- x$table
  out <- data.frame(Df = t$df, `Sum Sq` = t$ss, `Mean Sq` = t$ms, `F value` = t$F,
                    `Pr(>F)` = t$p, check.names = FALSE)
  out <- rbind(out, data.frame(Df = x$residual[["df"]], `Sum Sq` = x$residual[["ss"]],
                               `Mean Sq` = x$residual[["ms"]], `F value` = NA, `Pr(>F)` = NA,
                               check.names = FALSE))
  rownames(out) <- c(t$source, if (length(x$strata) > 1) "Error (residual)" else "Residuals")
  fmt <- format(out, digits = digits)
  fmt[is.na(out)] <- ""
  fmt[[" "]] <- c(.stars(t$p), "")
  print(fmt, right = TRUE)
  tested <- t[!t$is_error & t$error != "Residuals", ]
  if (nrow(tested)) {
    cat("\nTested against a higher error stratum:\n")
    for (i in seq_len(nrow(tested)))
      cat("  ", tested$source[i], " vs ", t$source[t$term == tested$error[i]], "\n", sep = "")
  }
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat(sprintf("\nGrand mean = %s   CV = %.2f %%   R-squared = %.4f\n",
              format(x$grand_mean, digits = digits + 2), x$cv, x$r_squared))
  invisible(x)
}

#' Least-squares means of a design
#'
#' Means of the levels of one or more factors, adjusted for the other terms
#' of the model: the other factors get equal weights and covariates are held
#' at their mean.
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param factors Name(s) of the factor(s) whose means are wanted.
#' @return A data frame with one row per level (or combination of levels): the
#'   least-squares mean, its standard error, the number of observations and
#'   the raw mean and standard deviation.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha",
#'                     block = "Block", design = "rcbd")
#' ls_means(fit, "Nitrogen_kg_ha")
#' @export
ls_means <- function(fit, factors) {
  if (!inherits(fit, "ad_anova")) stop("'fit' must come from design_anova()", call. = FALSE)
  d <- fit$data
  facs <- names(d)[vapply(d, is.factor, NA)]
  bad <- setdiff(factors, facs)
  if (length(bad)) stop("not a factor of the model: ", paste(bad, collapse = ", "), call. = FALSE)
  others <- setdiff(facs, factors)
  lv <- lapply(d[facs], levels)
  target <- expand.grid(lv[factors], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  rest <- if (length(others)) expand.grid(lv[others], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  qf <- fit$qr_fit
  piv <- qf$qr$pivot[seq_len(qf$rank)]
  beta <- qf$coefficients[piv]
  Rinv <- backsolve(qf$qr$qr[seq_len(qf$rank), seq_len(qf$rank), drop = FALSE], diag(qf$rank))
  XtXi <- Rinv %*% t(Rinv)
  mse <- fit$residual[["ms"]]
  tt <- stats::delete.response(fit$model_terms)
  out <- lapply(seq_len(nrow(target)), function(i) {
    nd <- if (is.null(rest)) target[i, , drop = FALSE]
          else cbind(rest, target[rep(i, nrow(rest)), , drop = FALSE], row.names = NULL)
    for (f in facs) nd[[f]] <- factor(nd[[f]], levels = lv[[f]])
    for (x in fit$roles$covariate) nd[[x]] <- 0
    Xn <- stats::model.matrix(tt, nd, contrasts.arg = .contrasts_for(facs))
    xbar <- colMeans(Xn)[piv]
    sel <- rep(TRUE, nrow(d))
    for (f in factors) sel <- sel & d[[f]] == target[i, f]
    obs <- fit$y[sel]
    data.frame(target[i, , drop = FALSE],
               mean = sum(xbar * beta),
               se = sqrt(max(0, drop(t(xbar) %*% XtXi %*% xbar)) * mse),
               n = length(obs),
               raw_mean = if (length(obs)) mean(obs) else NA_real_,
               sd = if (length(obs) > 1) stats::sd(obs) else NA_real_,
               stringsAsFactors = FALSE, row.names = NULL)
  })
  out <- do.call(rbind, out)
  # a term tested against a higher stratum: SE from that stratum (balanced formula)
  if (length(factors) == 1) {
    e <- fit$table$error[fit$table$term == factors]
    if (length(e) && e != "Residuals") out$se <- sqrt(fit$strata[[e]]$ms / out$n)
  }
  out
}
