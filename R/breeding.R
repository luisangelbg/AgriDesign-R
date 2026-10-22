# Plant-breeding analyses: genotype x environment interaction (AMMI and the
# stability statistics) and mating designs (line x tester, Griffing's
# diallels and the North Carolina designs I, II and III).
# Every analysis needs balanced data. Sums of squares come from sequential
# least-squares fits (QR), which give the classical closed forms for
# balanced data; the effects of the diallels come from a least-squares fit
# with Griffing's restrictions imposed exactly (null space of the
# restrictions), so their standard errors are exact.

# ---------- shared helpers ----------

# Factor with one level per combination of its arguments, in order of appearance.
.bp_fac <- function(...) {
  x <- paste(..., sep = "\r")
  factor(x, levels = unique(x))
}

.bp_ind <- function(f) {
  f <- droplevels(as.factor(f))
  outer(as.integer(f), seq_len(nlevels(f)), "==") + 0
}

# Sequential (type I) sums of squares of a list of terms (factors or matrices),
# entered after the intercept.
.bp_seq <- function(y, terms) {
  X <- matrix(1, length(y), 1)
  q <- qr(X); rk <- q$rank; rss <- sum(qr.resid(q, y)^2)
  df <- ss <- numeric(length(terms))
  for (k in seq_along(terms)) {
    M <- terms[[k]]
    if (!is.matrix(M)) M <- .bp_ind(M)
    X <- cbind(X, M)
    q <- qr(X); r2 <- sum(qr.resid(q, y)^2)
    df[k] <- q$rank - rk; ss[k] <- max(0, rss - r2)
    rk <- q$rank; rss <- r2
  }
  list(df = stats::setNames(df, names(terms)), ss = stats::setNames(ss, names(terms)),
       rss = rss, df_res = length(y) - rk)
}

.bp_rows <- function(source, df, ss) {
  data.frame(source = source, df = df, ss = ss, ms = ifelse(df > 0, ss / df, NA_real_),
             F = NA_real_, p = NA_real_, against = "", stringsAsFactors = FALSE)
}

# F test of the rows `src` against the mean square of row `den` of the same table.
.bp_against <- function(tab, src, den) {
  j <- match(den, tab$source)
  if (is.na(j)) return(tab)
  for (s in src) {
    i <- match(s, tab$source)
    if (is.na(i)) next
    ok <- isTRUE(tab$df[i] > 0 && tab$df[j] > 0 && tab$ms[j] > 0)
    if (!ok) next
    tab$F[i] <- tab$ms[i] / tab$ms[j]
    tab$p[i] <- stats::pf(tab$F[i], tab$df[i], tab$df[j], lower.tail = FALSE)
    tab$against[i] <- trimws(den)
  }
  tab
}

# Standard error of a linear combination of mean squares, Var(MS) = 2 MS^2 / df.
.bp_lin_se <- function(coef, ms, df) sqrt(sum(coef^2 * 2 * ms^2 / df))

.bp_signif <- function(est, se, df) {
  se <- rep_len(se, length(est))
  t <- ifelse(se > 0, est / se, NA_real_)
  p <- if (is.finite(df) && df > 0) 2 * stats::pt(-abs(t), df) else rep(NA_real_, length(t))
  list(t = t, p = p)
}

.bp_print_anova <- function(tab, digits) {
  out <- data.frame(Df = tab$df, `Sum Sq` = tab$ss, `Mean Sq` = tab$ms, `F value` = tab$F,
                    `Pr(>F)` = tab$p, check.names = FALSE)
  fmt <- format(out, digits = digits)
  fmt[is.na(out)] <- ""
  fmt[[" "]] <- .stars(tab$p)
  rownames(fmt) <- tab$source
  print(fmt, right = TRUE)
  x <- tab[tab$against != "" & !tab$against %in% c("Residuals", "Error", "Pooled error"), , drop = FALSE]
  if (nrow(x)) {
    cat("\nTested against another mean square:\n")
    for (i in seq_len(nrow(x))) cat("  ", trimws(x$source[i]), " vs ", x$against[i], "\n", sep = "")
  }
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

.bp_print_df <- function(df, digits) {
  num <- vapply(df, is.numeric, NA)
  out <- df
  out[num] <- lapply(df[num], function(v) ifelse(is.na(v), "", format(v, digits = digits)))
  print(out, row.names = FALSE, right = TRUE)
}

.bp_input <- function(data, response, cols) {
  if (!is.data.frame(data)) stop("'data' must be a data frame", call. = FALSE)
  if (!is.character(response) || length(response) != 1)
    stop("'response' must be the name of one column", call. = FALSE)
  for (a in cols) if (!is.character(a) || length(a) != 1)
    stop("column arguments must be single column names", call. = FALSE)
  .check_columns(data, c(response, cols))
  d <- as.data.frame(data)[c(response, cols)]
  if (!is.numeric(d[[response]])) stop("the response '", response, "' must be numeric", call. = FALSE)
  d
}

.bp_check_inbreeding <- function(f) {
  if (!is.numeric(f) || length(f) != 1 || is.na(f) || f < 0 || f > 1)
    stop("'inbreeding' must be a number between 0 and 1", call. = FALSE)
}

# ---------- genotype x environment table ----------

# Two-way table of genotype means from plot data (or from one mean per cell),
# with the pooled error of the trial on the plot scale.
.gxe_input <- function(data, response, genotype, environment, block, mse, df_error, reps) {
  d <- .bp_input(data, response, c(genotype, environment, block))
  ok <- stats::complete.cases(d)
  d <- d[ok, , drop = FALSE]
  G <- factor(as.character(d[[genotype]]), levels = .levels_of(d[[genotype]]))
  E <- factor(as.character(d[[environment]]), levels = .levels_of(d[[environment]]))
  g <- nlevels(G); e <- nlevels(E)
  if (g < 3) stop("at least three genotypes are needed", call. = FALSE)
  if (e < 3) stop("at least three environments are needed", call. = FALSE)
  n <- table(G, E)
  if (any(n == 0)) {
    w <- which(n == 0, arr.ind = TRUE)
    lab <- paste(rownames(n)[w[, 1]], colnames(n)[w[, 2]], sep = " x ")
    stop("every genotype must be present in every environment; missing: ",
         paste(utils::head(lab, 5), collapse = ", "), if (length(lab) > 5) ", ...", call. = FALSE)
  }
  if (length(unique(as.vector(n))) != 1)
    stop("the number of replicates must be the same in every genotype x environment cell", call. = FALSE)
  k <- n[1, 1]
  y <- d[[response]]
  Y <- matrix(tapply(y, list(G, E), mean), g, e, dimnames = list(levels(G), levels(E)))
  out <- list(Y = Y, n = length(y), dropped = sum(!ok), replicated = k > 1)
  A <- .gxe_additive(Y)
  if (k > 1) {
    if (!is.null(mse)) stop("'mse' is only used when the data hold one value (a mean) per ",
                            "genotype x environment cell", call. = FALSE)
    terms <- list(ENV = E)
    if (!is.null(block)) {
      B <- .bp_fac(E, d[[block]])
      if (any(table(B, G) != 1))
        stop("every block must hold each genotype exactly once (complete blocks within ",
             "each environment)", call. = FALSE)
      terms$REP <- B
    }
    terms$GEN <- G
    terms$GxE <- .bp_fac(G, E)
    sq <- .bp_seq(y, terms)
    out$r <- k; out$dfe <- sq$df_res; out$mse <- sq$rss / sq$df_res
    out$base <- .bp_rows(c("ENV", if (!is.null(block)) "REP(ENV)", "GEN", "GxE"),
                         unname(sq$df), unname(sq$ss))
    out$total <- c(df = length(y) - 1, ss = sum((y - mean(y))^2))
  } else {
    if (!is.null(block)) stop("with one value per genotype x environment cell there are no ",
                              "replicates: drop 'block' and give 'mse' and 'df_error'", call. = FALSE)
    if (!is.numeric(reps) || length(reps) != 1 || reps < 1) stop("'reps' must be a positive number", call. = FALSE)
    out$r <- reps
    if (!is.null(mse)) {
      if (is.null(df_error)) stop("'df_error' is needed together with 'mse'", call. = FALSE)
      if (!is.numeric(mse) || mse <= 0 || !is.numeric(df_error) || df_error < 1)
        stop("'mse' must be positive and 'df_error' at least 1", call. = FALSE)
      out$mse <- mse; out$dfe <- df_error
    } else { out$mse <- NA_real_; out$dfe <- NA_real_ }
    out$base <- .bp_rows(c("ENV", "GEN", "GxE"), c(e - 1, g - 1, (g - 1) * (e - 1)),
                         reps * c(g * sum((A$em - A$grand)^2), e * sum((A$gm - A$grand)^2), sum(A$Z^2)))
    out$total <- NULL
  }
  out
}

.gxe_additive <- function(Y) {
  grand <- mean(Y); gm <- rowMeans(Y); em <- colMeans(Y)
  list(grand = grand, gm = gm, em = em, Z = Y - outer(gm, em, "+") + grand)
}

# ---------- AMMI ----------

#' AMMI analysis of a multi-environment trial
#'
#' Additive main effects and multiplicative interaction model (Gollob 1968;
#' Gauch 1988): the analysis of variance of genotypes and environments, and
#' the singular value decomposition of the genotype x environment
#' interaction into principal axes (IPCA), each tested with the degrees of
#' freedom of Gollob (1968).
#'
#' The data must be balanced: every genotype in every environment with the
#' same number of replicates. With plot data the error comes from the trial
#' itself (with `block`, blocks within environments are taken out and
#' environments are tested against them). With one mean per cell, give the
#' pooled error of the trial (`mse`, on the plot scale, with `df_error`) and
#' the number of replicates behind each mean (`reps`).
#'
#' Interaction scores are scaled by the square root of the singular values
#' (symmetric scaling), so that the product of a genotype score and an
#' environment score is their interaction. The sign of every axis is fixed so
#' that the environment with the largest absolute score is positive.
#' The AMMI stability value (ASV) of Purchase et al. (2000) is
#' `sqrt((SS1 / SS2 * IPCA1)^2 + IPCA2^2)`; the yield stability index (YSI)
#' adds its rank to the rank of the mean (highest mean = 1).
#'
#' @param data A data frame in long format.
#' @param response Name of the numeric response column.
#' @param genotype,environment Names of the genotype and environment columns.
#' @param block Optional name of the block (replicate) column; blocks are
#'   taken as nested within environments.
#' @param mse,df_error Error mean square (plot scale) and its degrees of
#'   freedom, only when the data hold one mean per cell.
#' @param reps Number of replicates behind each mean, when the data hold one
#'   mean per cell.
#' @param n_axes Number of interaction axes shown in the analysis of
#'   variance (the rest is pooled as the residual interaction); by default
#'   all of them.
#' @return An object of class `"ad_ammi"`: the analysis of variance
#'   (`anova`), the interaction axes (`axes`: singular value, sum of
#'   squares, percentage of the interaction, Gollob degrees of freedom and F
#'   test), genotype and environment means with their scores (`genotypes`,
#'   `environments`, the first with ASV and YSI), the table of means and the
#'   interaction residuals.
#' @references Gollob, H. F. (1968) Psychometrika 33, 73-115.
#'   Gauch, H. G. (1988) Biometrics 44, 705-715.
#'   Purchase, J. L., Hatting, H. and van Deventer, C. S. (2000) South
#'   African Journal of Plant and Soil 17, 101-107.
#' @examples
#' fit <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
#' fit
#' plot(fit)
#' plot(fit, type = "ammi2")
#' @export
ammi <- function(data, response, genotype, environment, block = NULL, mse = NULL,
                 df_error = NULL, reps = 1, n_axes = NULL) {
  tw <- .gxe_input(data, response, genotype, environment, block, mse, df_error, reps)
  Y <- tw$Y; g <- nrow(Y); e <- ncol(Y); r <- tw$r
  A <- .gxe_additive(Y)
  sv <- svd(A$Z)
  m <- min(g - 1, e - 1)
  d <- sv$d[seq_len(m)]
  U <- sv$u[, seq_len(m), drop = FALSE]; V <- sv$v[, seq_len(m), drop = FALSE]
  for (k in seq_len(m)) {
    j <- which.max(abs(V[, k]))
    if (V[j, k] < 0) { U[, k] <- -U[, k]; V[, k] <- -V[, k] }
  }
  if (is.null(n_axes)) n_axes <- m
  if (!is.numeric(n_axes) || length(n_axes) != 1 || n_axes < 1 || n_axes > m || n_axes != round(n_axes))
    stop("'n_axes' must be a whole number between 1 and ", m, call. = FALSE)
  ss_ge <- r * sum(A$Z^2)
  ss_ax <- r * d^2
  df_ax <- g + e - 1 - 2 * seq_len(m)
  axes <- data.frame(axis = paste0("IPCA", seq_len(m)), singular_value = d, ss = ss_ax,
                     pct = 100 * ss_ax / ss_ge, cum_pct = 100 * cumsum(ss_ax) / ss_ge,
                     df = df_ax, ms = ss_ax / df_ax, F = NA_real_, p = NA_real_, stringsAsFactors = FALSE)
  has_err <- is.finite(tw$mse)
  if (has_err) {
    axes$F <- axes$ms / tw$mse
    axes$p <- stats::pf(axes$F, axes$df, tw$dfe, lower.tail = FALSE)
  }
  k <- seq_len(n_axes)
  tab <- rbind(tw$base, .bp_rows(paste0("  IPCA", k), df_ax[k], ss_ax[k]))
  if (n_axes < m) tab <- rbind(tab, .bp_rows("  Residual GxE", sum(df_ax[-k]), sum(ss_ax[-k])))
  if (has_err) tab <- rbind(tab, .bp_rows("Residuals", tw$dfe, tw$mse * tw$dfe))
  if (!is.null(tw$total)) tab <- rbind(tab, .bp_rows("Total", tw$total[["df"]], tw$total[["ss"]]))
  if (has_err) {
    tab <- .bp_against(tab, "ENV", if ("REP(ENV)" %in% tab$source) "REP(ENV)" else "Residuals")
    tab <- .bp_against(tab, c("REP(ENV)", "GEN", "GxE", paste0("  IPCA", k), "  Residual GxE"), "Residuals")
  }
  sc <- function(M) M %*% diag(sqrt(d), m, m)
  Gs <- sc(U); Es <- sc(V)
  colnames(Gs) <- colnames(Es) <- axes$axis
  gen <- data.frame(genotype = rownames(Y), mean = unname(A$gm), Gs, stringsAsFactors = FALSE,
                    row.names = NULL, check.names = FALSE)
  env <- data.frame(environment = colnames(Y), mean = unname(A$em), Es, stringsAsFactors = FALSE,
                    row.names = NULL, check.names = FALSE)
  if (m >= 2) {
    gen$asv <- sqrt((ss_ax[1] / ss_ax[2] * Gs[, 1])^2 + Gs[, 2]^2)
    gen$rank_asv <- rank(signif(gen$asv, 12), ties.method = "min")
    gen$rank_mean <- rank(-signif(gen$mean, 12), ties.method = "min")
    gen$ysi <- gen$rank_asv + gen$rank_mean
  }
  structure(list(anova = tab, axes = axes, genotypes = gen, environments = env, means = Y,
                 interaction = A$Z, grand_mean = A$grand, r = r, mse = tw$mse, df_error = tw$dfe,
                 n_axes = n_axes, response = response, n = tw$n, dropped = tw$dropped,
                 replicated = tw$replicated),
            class = "ad_ammi")
}

#' @export
print.ad_ammi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nAMMI analysis (additive main effects and multiplicative interaction)\n")
  cat("Response: ", x$response, "   ", nrow(x$means), " genotypes x ", ncol(x$means), " environments, ",
      if (x$replicated) paste(x$r, "replicate(s) per cell") else paste("means of", x$r, "replicate(s)"),
      if (x$dropped) paste0(", ", x$dropped, " incomplete row(s) dropped"), "\n\n", sep = "")
  .bp_print_anova(x$anova, digits)
  cat("\nInteraction principal components (degrees of freedom of Gollob 1968):\n")
  ax <- x$axes[, c("axis", "singular_value", "pct", "cum_pct", "df")]
  names(ax) <- c("Axis", "Singular value", "% of GxE", "Cumulative %", "Df")
  .bp_print_df(ax, digits)
  cat("\nGenotypes (scores scaled by the square root of the singular values):\n")
  .bp_print_df(x$genotypes, digits)
  cat("\nEnvironments:\n")
  .bp_print_df(x$environments, digits)
  if (!is.null(x$genotypes$asv))
    cat("\nASV: AMMI stability value (Purchase et al. 2000); YSI = rank of ASV + rank of the mean.\n")
  invisible(x)
}

#' Biplots of an AMMI analysis
#'
#' The AMMI1 biplot puts the means of genotypes and environments on the
#' horizontal axis and their first interaction score on the vertical axis;
#' the AMMI2 biplot shows the first two interaction axes, with the
#' environments as vectors from the origin.
#'
#' @param x An `"ad_ammi"` object from [ammi()].
#' @param type `"ammi1"` (mean vs IPCA1) or `"ammi2"` (IPCA1 vs IPCA2).
#' @param col Two colours: genotypes and environments.
#' @param labels Write the names of genotypes and environments.
#' @param main,xlab,ylab Titles; `NULL` builds them from the analysis.
#' @param cex Size of the labels.
#' @param ... Further graphical parameters passed to [graphics::plot.window()].
#' @return Invisibly, a data frame with the coordinates of every point.
#' @examples
#' fit <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
#' plot(fit, type = "ammi2")
#' @export
plot.ad_ammi <- function(x, type = c("ammi1", "ammi2"), col = NULL, labels = TRUE,
                         main = NULL, xlab = NULL, ylab = NULL, cex = 0.8, ...) {
  type <- match.arg(tolower(type[1]), c("ammi1", "ammi2"))
  col <- if (is.null(col)) agri_colors(2) else rep_len(col, 2)
  G <- x$genotypes; E <- x$environments; ax <- x$axes
  lab_ax <- function(k) sprintf("IPCA%d (%.1f %%)", k, ax$pct[k])
  if (type == "ammi1") {
    xg <- G$mean; yg <- G$IPCA1; xe <- E$mean; ye <- E$IPCA1
    xl <- paste("Mean of", x$response); yl <- lab_ax(1)
  } else {
    if (nrow(ax) < 2) stop("the AMMI2 biplot needs two interaction axes", call. = FALSE)
    xg <- G$IPCA1; yg <- G$IPCA2; xe <- E$IPCA1; ye <- E$IPCA2
    xl <- lab_ax(1); yl <- lab_ax(2)
  }
  pad <- function(v) { r <- range(v); w <- diff(r); if (w == 0) w <- 1; r + c(-1, 1) * w * 0.12 }
  xlim <- pad(c(xg, xe, if (type == "ammi2") 0)); ylim <- pad(c(yg, ye, 0))
  op <- .style_par(); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim, asp = if (type == "ammi2") 1 else NA, ...)
  graphics::abline(h = 0, v = if (type == "ammi1") x$grand_mean else 0, lty = 2, col = "#9aa0a6")
  if (type == "ammi2")
    graphics::arrows(0, 0, xe, ye, length = 0.07, angle = 22, col = .lighten(col[2], 0.25), lwd = 1.3)
  graphics::points(xg, yg, pch = 21, bg = col[1], col = .darken(col[1]), cex = 1.2)
  graphics::points(xe, ye, pch = 24, bg = col[2], col = .darken(col[2]), cex = 1.2)
  if (labels) {
    graphics::text(xg, yg, G$genotype, pos = 3, cex = cex, col = .darken(col[1], 0.35))
    graphics::text(xe, ye, E$environment, pos = 1, cex = cex, col = .darken(col[2], 0.35), font = 2)
  }
  graphics::axis(1); graphics::axis(2); graphics::box()
  graphics::title(main = if (is.null(main)) paste(if (type == "ammi1") "AMMI1" else "AMMI2", "biplot") else main,
                  xlab = if (is.null(xlab)) xl else xlab, ylab = if (is.null(ylab)) yl else ylab,
                  font.main = 1, cex.main = 1.05)
  graphics::legend("topleft", c("Genotypes", "Environments"), pch = c(21, 24), pt.bg = col,
                   col = .darken(col), bty = "n", cex = 0.8, inset = 0.01)
  .note(sprintf("IPCA1 + IPCA2 = %.1f %% of GxE", min(100, sum(ax$pct[seq_len(min(2, nrow(ax)))]))))
  invisible(data.frame(label = c(G$genotype, E$environment),
                       type = rep(c("genotype", "environment"), c(length(xg), length(xe))),
                       x = c(xg, xe), y = c(yg, ye), stringsAsFactors = FALSE))
}

# ---------- stability ----------

#' Stability statistics of a multi-environment trial
#'
#' Parametric and rank-based stability statistics computed from the table of
#' genotype means in every environment and the pooled error of the trial:
#' \describe{
#'   \item{`"eberhart"`}{Eberhart and Russell (1966): regression coefficient
#'     `b` on the environmental index (with its standard error from the
#'     genotype's own deviations, and a t test of `b = 1`), deviation mean
#'     square `s2d` (mean scale, minus the error of a mean) with its F test
#'     against the pooled error, and their analysis of variance.}
#'   \item{`"finlay"`}{Finlay and Wilkinson (1963): regression on the
#'     environment mean, optionally on the log scale.}
#'   \item{`"shukla"`}{Shukla's (1972) stability variance `sigma2` and the
#'     same after removing the linear regression on the environmental index
#'     (`s2`), with F tests against the pooled error (plot scale).}
#'   \item{`"wricke"`}{Wricke's (1962) ecovalence (plot scale) and its
#'     share of the interaction.}
#'   \item{`"kang"`}{Kang's (1993) yield-stability statistic YS: the rank of
#'     the mean (highest = g), adjusted by the LSD, plus the stability
#'     rating from Shukla's test (-8, -4, -2 at 0.01, 0.05, 0.10).}
#'   \item{`"huehn"`}{Huehn's (1979) rank statistics S1 and S2 on the ranks
#'     corrected for the genotype effects, with the Z tests of Nassar and
#'     Huehn (1987).}
#' }
#' The common analysis of variance follows Shukla (1972): the interaction is
#' split into the heterogeneity of the regressions and the residual, with
#' (g - 1)(e - 2) degrees of freedom (Eberhart and Russell printed g(e - 2),
#' but e - 2 of them belong to environments and are identically zero).
#'
#' @inheritParams ammi
#' @param method One or more of `"eberhart"`, `"finlay"`, `"shukla"`,
#'   `"wricke"`, `"kang"` and `"huehn"` (default: all).
#' @param alpha Significance level of the LSD in Kang's statistic and of the
#'   critical values of Huehn's tests.
#' @param log_finlay Fit the regression of Finlay and Wilkinson on the
#'   natural logarithm of the means, as in their paper.
#' @return An object of class `"ad_stability"`: the analysis of variance
#'   (`anova`) and one element per method.
#' @references Eberhart, S. A. and Russell, W. A. (1966) Crop Science 6, 36-40.
#'   Finlay, K. W. and Wilkinson, G. N. (1963) Australian Journal of
#'   Agricultural Research 14, 742-754.
#'   Shukla, G. K. (1972) Heredity 29, 237-245.
#'   Wricke, G. (1962) Zeitschrift fur Pflanzenzuchtung 47, 92-96.
#'   Kang, M. S. (1993) Agronomy Journal 85, 754-757.
#'   Huehn, M. (1979) EDV in Medizin und Biologie 10, 112-117.
#'   Nassar, R. and Huehn, M. (1987) Biometrics 43, 45-53.
#' @examples
#' st <- stability(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
#' st
#' st$kang
#' @export
stability <- function(data, response, genotype, environment, block = NULL, mse = NULL,
                      df_error = NULL, reps = 1,
                      method = c("eberhart", "finlay", "shukla", "wricke", "kang", "huehn"),
                      alpha = 0.05, log_finlay = FALSE) {
  method <- match.arg(method, several.ok = TRUE)
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be between 0 and 1", call. = FALSE)
  tw <- .gxe_input(data, response, genotype, environment, block, mse, df_error, reps)
  Y <- tw$Y; g <- nrow(Y); e <- ncol(Y); r <- tw$r; mse <- tw$mse; dfe <- tw$dfe
  A <- .gxe_additive(Y)
  I <- A$em - A$grand
  sI2 <- sum(I^2)
  if (sI2 <= 0) stop("the environments have the same mean: the environmental index is zero", call. = FALSE)
  b <- drop(Y %*% I) / sI2
  dev <- Y - A$gm - outer(b, I)
  ssdev <- rowSums(dev^2)
  msdev <- ssdev / (e - 2)
  has_err <- is.finite(mse)
  mse_mean <- mse / r
  # analysis of variance (Shukla 1972), plot scale
  ss_het <- r * (sum(b^2) - g) * sI2
  ss_ge <- r * sum(A$Z^2)
  tab <- .bp_rows(c("Genotypes", "Environments", "GxE", "  Heterogeneity", "  Residual"),
                  c(g - 1, e - 1, (g - 1) * (e - 1), g - 1, (g - 1) * (e - 2)),
                  c(r * e * sum((A$gm - A$grand)^2), r * g * sI2, ss_ge, ss_het, ss_ge - ss_het))
  if (has_err) tab <- rbind(tab, .bp_rows("Pooled error", dfe, mse * dfe))
  tab <- rbind(tab, .bp_rows("Total", g * e - 1, r * sum((Y - A$grand)^2)))
  tab <- .bp_against(tab, "Genotypes", "GxE")
  tab <- .bp_against(tab, "  Heterogeneity", "  Residual")
  if (has_err) tab <- .bp_against(tab, c("Environments", "GxE", "  Residual"), "Pooled error")
  out <- list(anova = tab, method = method, means = Y, r = r, mse = mse, df_error = dfe,
              response = response, n = tw$n, dropped = tw$dropped, replicated = tw$replicated)
  gen <- rownames(Y)
  pdev <- if (has_err) stats::pf(msdev / mse_mean, e - 2, dfe, lower.tail = FALSE) else NA_real_

  if ("eberhart" %in% method) {
    se_b <- sqrt(msdev / sI2)
    tt <- .bp_signif(b - 1, se_b, e - 2)
    ss_tot <- rowSums((Y - A$gm)^2)
    er <- data.frame(genotype = gen, mean = unname(A$gm), b = b, se_b = se_b, t_b1 = tt$t, p_b1 = tt$p,
                     s2d = if (has_err) msdev - mse_mean else NA_real_,
                     F_dev = if (has_err) msdev / mse_mean else NA_real_, p_dev = pdev,
                     r2 = ifelse(ss_tot > 0, 1 - ssdev / ss_tot, NA_real_),
                     stringsAsFactors = FALSE, row.names = NULL)
    pooled <- sum(ssdev)
    et <- .bp_rows(c("Genotypes", "Env + (G x E)", "  Environments (linear)", "  G x E (linear)",
                     "  Pooled deviation"),
                   c(g - 1, g * (e - 1), 1, g - 1, (g - 1) * (e - 2)),
                   c(r * e * sum((A$gm - A$grand)^2), r * sum((Y - A$gm)^2), r * g * sI2, ss_het, r * pooled))
    if (has_err) et <- rbind(et, .bp_rows("Pooled error", dfe, mse * dfe))
    et <- .bp_against(et, c("Genotypes", "  G x E (linear)"), "  Pooled deviation")
    if (has_err) et <- .bp_against(et, "  Pooled deviation", "Pooled error")
    out$eberhart <- list(table = er, anova = et, index = stats::setNames(I, colnames(Y)))
  }
  if ("finlay" %in% method) {
    L <- Y
    if (log_finlay) {
      if (any(Y <= 0)) stop("the log scale of Finlay and Wilkinson needs positive means", call. = FALSE)
      L <- log(Y)
    }
    x <- colMeans(L); xc <- x - mean(x); sx <- sum(xc^2)
    bf <- drop((L - rowMeans(L)) %*% xc) / sx
    af <- rowMeans(L) - bf * mean(x)
    rf <- L - af - outer(bf, x)
    sf <- rowSums(rf^2)
    st <- rowSums((L - rowMeans(L))^2)
    out$finlay <- data.frame(genotype = gen, mean = unname(A$gm), intercept = unname(af), b = unname(bf),
                             se_b = sqrt(sf / (e - 2) / sx), r2 = ifelse(st > 0, 1 - sf / st, NA_real_),
                             stringsAsFactors = FALSE, row.names = NULL)
    attr(out$finlay, "log") <- log_finlay
  }
  W <- rowSums(A$Z^2)
  sig <- r * (g * (g - 1) * W - sum(W)) / ((g - 1) * (g - 2) * (e - 1))
  F_sig <- if (has_err) sig / mse else rep(NA_real_, g)
  p_sig <- if (has_err) stats::pf(F_sig, e - 1, dfe, lower.tail = FALSE) else rep(NA_real_, g)
  if ("shukla" %in% method) {
    s2 <- r * g / ((g - 2) * (e - 2)) * (ssdev - sum(ssdev) / (g * (g - 1)))
    out$shukla <- data.frame(genotype = gen, mean = unname(A$gm), sigma2 = unname(sig), F = unname(F_sig),
                             p = unname(p_sig), s2 = unname(s2),
                             F_s2 = if (has_err) unname(s2 / mse) else NA_real_,
                             p_s2 = if (has_err) stats::pf(unname(s2 / mse), e - 2, dfe, lower.tail = FALSE) else NA_real_,
                             stringsAsFactors = FALSE, row.names = NULL)
  }
  if ("wricke" %in% method)
    out$wricke <- data.frame(genotype = gen, mean = unname(A$gm), W = unname(r * W),
                             pct = unname(100 * W / sum(W)), stringsAsFactors = FALSE, row.names = NULL)
  if ("kang" %in% method) {
    if (has_err) {
      lsd <- stats::qt(1 - alpha / 2, dfe) * sqrt(2 * mse / (r * e))
      rk <- rank(signif(A$gm, 12), ties.method = "min")
      dm <- A$gm - mean(A$gm)
      adj <- sign(dm)
      adj[dm >= lsd] <- 2; adj[dm >= 2 * lsd] <- 3
      adj[dm <= -lsd] <- -2; adj[dm <= -2 * lsd] <- -3
      rating <- ifelse(p_sig < 0.01, -8, ifelse(p_sig < 0.05, -4, ifelse(p_sig < 0.10, -2, 0)))
      ys <- rk + adj + rating
      kt <- data.frame(genotype = gen, mean = unname(A$gm), rank = unname(rk), adjustment = unname(adj),
                       adjusted_rank = unname(rk + adj), sigma2 = unname(sig), rating = unname(rating),
                       ys = unname(ys), selected = unname(ys > mean(ys)), stringsAsFactors = FALSE,
                       row.names = NULL)
      out$kang <- list(table = kt, lsd = lsd, mean_ys = mean(ys), alpha = alpha)
    } else out$kang <- NULL
  }
  if ("huehn" %in% method) {
    Xc <- Y - A$gm + A$grand
    R <- apply(signif(Xc, 12), 2, rank)
    K <- g; N <- e
    S1 <- apply(R, 1, function(v) sum(abs(outer(v, v, "-")))) / (N * (N - 1))
    S2 <- apply(R, 1, stats::var)
    E1 <- (K^2 - 1) / (3 * K)
    V1 <- (K^2 - 1) * ((K^2 - 4) * (N + 3) + 30) / (45 * K^2 * N * (N - 1))
    E2 <- (K^2 - 1) / 12
    V2 <- (K^2 - 1) * (2 * (K^2 - 4) * (N - 1) + 5 * (K^2 - 1)) / (360 * N * (N - 1))
    Z1 <- (S1 - E1)^2 / V1; Z2 <- (S2 - E2)^2 / V2
    ht <- data.frame(genotype = gen, mean = unname(A$gm), rank_mean = rank(signif(A$gm, 12)),
                     S1 = unname(S1), Z1 = unname(Z1), p1 = stats::pchisq(unname(Z1), 1, lower.tail = FALSE),
                     S2 = unname(S2), Z2 = unname(Z2), p2 = stats::pchisq(unname(Z2), 1, lower.tail = FALSE),
                     stringsAsFactors = FALSE, row.names = NULL)
    out$huehn <- list(table = ht, ranks = R,
                      expectation = c(E1 = E1, V1 = V1, E2 = E2, V2 = V2),
                      sum_z = c(Z1 = sum(Z1), Z2 = sum(Z2)),
                      p_sum = c(Z1 = stats::pchisq(sum(Z1), K, lower.tail = FALSE),
                                Z2 = stats::pchisq(sum(Z2), K, lower.tail = FALSE)),
                      critical = c(individual = stats::qchisq(1 - alpha / K, 1),
                                   sum = stats::qchisq(1 - alpha, K)))
  }
  structure(out, class = "ad_stability")
}

#' @export
print.ad_stability <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nStability analysis of a multi-environment trial\n")
  cat("Response: ", x$response, "   ", nrow(x$means), " genotypes x ", ncol(x$means), " environments, ",
      if (x$replicated) paste(x$r, "replicate(s) per cell") else paste("means of", x$r, "replicate(s)"), "\n", sep = "")
  if (!is.finite(x$mse)) cat("No error mean square: the tests against the pooled error are not available.\n")
  cat("\nAnalysis of variance (plot scale; Shukla 1972):\n")
  .bp_print_anova(x$anova, digits)
  if (!is.null(x$eberhart)) {
    cat("\nEberhart and Russell (1966): regression on the environmental index\n")
    .bp_print_anova(x$eberhart$anova, digits)
    cat("\n")
    .bp_print_df(x$eberhart$table, digits)
    cat("b: SE from the genotype's own deviations (e - 2 df), t tests b = 1;",
        "s2d = MS deviations - error of a mean.\n")
  }
  if (!is.null(x$finlay)) {
    cat("\nFinlay and Wilkinson (1963): regression on the environment mean",
        if (isTRUE(attr(x$finlay, "log"))) "(log scale)" else "", "\n")
    .bp_print_df(x$finlay, digits)
  }
  if (!is.null(x$shukla)) {
    cat("\nShukla (1972): stability variance sigma2 and s2 (heterogeneity removed), plot scale\n")
    .bp_print_df(x$shukla, digits)
  }
  if (!is.null(x$wricke)) {
    cat("\nWricke (1962): ecovalence W (plot scale) and its share of the interaction (%)\n")
    .bp_print_df(x$wricke, digits)
  }
  if (!is.null(x$kang)) {
    cat("\nKang (1993): yield-stability statistic YS\n")
    .bp_print_df(x$kang$table, digits)
    cat(sprintf("LSD (alpha = %g) = %s   mean YS = %s; selected: YS above the mean\n", x$kang$alpha,
                format(x$kang$lsd, digits = digits), format(x$kang$mean_ys, digits = digits)))
  } else if ("kang" %in% x$method) cat("\nKang (1993): needs the error mean square of the trial.\n")
  if (!is.null(x$huehn)) {
    h <- x$huehn
    cat("\nHuehn (1979): rank statistics with the Z tests of Nassar and Huehn (1987)\n")
    .bp_print_df(h$table, digits)
    cat(sprintf("E(S1) = %s  Var(S1) = %s  E(S2) = %s  Var(S2) = %s\n", format(h$expectation[1], digits = digits),
                format(h$expectation[2], digits = digits), format(h$expectation[3], digits = digits),
                format(h$expectation[4], digits = digits)))
    cat(sprintf("Sum of Z1 = %s (p = %s)   Sum of Z2 = %s (p = %s); chi-square critical values: %s (each Z, Bonferroni), %s (sum)\n",
                format(h$sum_z[1], digits = digits), format.pval(h$p_sum[1], digits = 3),
                format(h$sum_z[2], digits = digits), format.pval(h$p_sum[2], digits = 3),
                format(h$critical[1], digits = digits), format(h$critical[2], digits = digits)))
  }
  invisible(x)
}

# ---------- line x tester ----------

#' Line x tester analysis
#'
#' Analysis of a line x tester mating design (Kempthorne 1957; Singh and
#' Chaudhary 1985): the analysis of variance of the entries in randomised
#' complete blocks, with the crosses partitioned into lines, testers and
#' their interaction; general combining ability of lines and testers,
#' specific combining ability of the crosses, their standard errors, the
#' proportional contribution of every source to the variation among crosses
#' and the genetic components.
#'
#' A row is a cross when it has both a line and a tester, a line parent when
#' the tester is missing (`NA` or empty) and a tester parent when the line is
#' missing; parents are optional. Rows with neither (for instance checks) are
#' not used.
#'
#' Lines and testers are tested against the error and, in the random model
#' of Kempthorne, against the line x tester interaction (`anova_lxt`). The
#' average covariance of half sibs is
#' `[((l-1) MSl + (t-1) MSt) / (l+t-2) - MSlxt] (l+t-2) / (r (2lt - l - t))`,
#' the moment estimator under equal variances of lines and testers;
#' `sigma2_A = 4 CovHS / (1 + F)` and `sigma2_D = 4 sigma2_SCA / (1 + F)^2`
#' are given for F = 0 and F = 1. The standard errors of the effects are
#' the usual ones, `sqrt(MSE / (r t))` for lines, `sqrt(MSE / (r l))` for
#' testers and `sqrt(MSE / r)` for crosses.
#'
#' @param data A data frame in long format.
#' @param response Name of the numeric response column.
#' @param line,tester Names of the line (female) and tester (male) columns.
#' @param block Name of the block column (randomised complete blocks); when
#'   `NULL` the entries are taken as completely randomised.
#' @return An object of class `"ad_line_tester"`.
#' @references Kempthorne, O. (1957) An Introduction to Genetic Statistics.
#'   Wiley. Singh, R. K. and Chaudhary, B. D. (1985) Biometrical Methods in
#'   Quantitative Genetic Analysis. Kalyani.
#' @examples
#' lt <- line_tester(maize_lxt, "Yield_t_ha", "Line", "Tester", block = "Block")
#' lt
#' @export
line_tester <- function(data, response, line, tester, block = NULL) {
  d <- .bp_input(data, response, c(line, tester, block))
  ok <- !is.na(d[[response]]) & (if (is.null(block)) TRUE else !is.na(d[[block]]))
  d <- d[ok, , drop = FALSE]
  lv <- as.character(d[[line]]); tv <- as.character(d[[tester]])
  lv[!is.na(lv) & trimws(lv) == ""] <- NA; tv[!is.na(tv) & trimws(tv) == ""] <- NA
  checks <- is.na(lv) & is.na(tv)
  d <- d[!checks, , drop = FALSE]; lv <- lv[!checks]; tv <- tv[!checks]
  kind <- ifelse(!is.na(lv) & !is.na(tv), "cross", ifelse(is.na(tv), "line", "tester"))
  cr <- kind == "cross"
  if (!any(cr)) stop("no crosses: a cross needs both a line and a tester", call. = FALSE)
  L <- .levels_of(d[[line]][cr]); Tt <- .levels_of(d[[tester]][cr])
  l <- length(L); t <- length(Tt)
  if (l < 2 || t < 2) stop("at least two lines and two testers are needed", call. = FALSE)
  tab_c <- table(factor(lv[cr], levels = L), factor(tv[cr], levels = Tt))
  if (any(tab_c == 0)) {
    w <- which(tab_c == 0, arr.ind = TRUE)
    stop("every line must be crossed with every tester; missing: ",
         paste(utils::head(paste(L[w[, 1]], Tt[w[, 2]], sep = " x "), 5), collapse = ", "), call. = FALSE)
  }
  entry <- ifelse(kind == "cross", paste(lv, tv, sep = " x "), ifelse(kind == "line", lv, tv))
  key <- paste(kind, entry, sep = "\r")
  Ent <- factor(key, levels = unique(key))
  y <- d[[response]]
  if (!is.null(block)) {
    B <- factor(as.character(d[[block]]))
    if (nlevels(B) < 2) stop("at least two blocks are needed", call. = FALSE)
    if (any(table(B, Ent) != 1)) stop("every entry must appear exactly once in every block", call. = FALSE)
    r <- nlevels(B)
    sq <- .bp_seq(y, list(Blocks = B, Entries = Ent))
  } else {
    cnt <- table(Ent)
    if (length(unique(cnt)) != 1 || cnt[1] < 2)
      stop("every entry needs the same number of replicates (at least two)", call. = FALSE)
    r <- unname(cnt[1])
    sq <- .bp_seq(y, list(Entries = Ent))
  }
  dfe <- sq$df_res; mse <- sq$rss / dfe
  if (dfe < 1) stop("no degrees of freedom left for the error", call. = FALSE)
  em <- tapply(y, Ent, mean)
  ekind <- sub("\r.*$", "", names(em))
  M <- matrix(tapply(y[cr], list(factor(lv[cr], levels = L), factor(tv[cr], levels = Tt)), mean), l, t,
              dimnames = list(L, Tt))
  gm <- mean(M); lm_ <- rowMeans(M); tm <- colMeans(M)
  ss_cross <- r * sum((M - gm)^2)
  ss_l <- r * t * sum((lm_ - gm)^2); ss_t <- r * l * sum((tm - gm)^2)
  ss_lt <- ss_cross - ss_l - ss_t
  pm <- em[ekind != "cross"]; pk <- ekind[ekind != "cross"]
  np <- length(pm)
  rows <- list()
  if (!is.null(block)) rows$b <- .bp_rows("Blocks", sq$df[["Blocks"]], sq$ss[["Blocks"]])
  rows$e <- .bp_rows("Entries", sq$df[["Entries"]], sq$ss[["Entries"]])
  if (np >= 2) {
    ss_p <- r * sum((pm - mean(pm))^2)
    rows$p <- .bp_rows("  Parents", np - 1, ss_p)
    nl <- sum(pk == "line"); nt <- sum(pk == "tester")
    if (nl >= 1 && nt >= 1) {
      sl <- if (nl > 1) r * sum((pm[pk == "line"] - mean(pm[pk == "line"]))^2) else 0
      st <- if (nt > 1) r * sum((pm[pk == "tester"] - mean(pm[pk == "tester"]))^2) else 0
      if (nl > 1) rows$pl <- .bp_rows("    Lines (parents)", nl - 1, sl)
      if (nt > 1) rows$pt <- .bp_rows("    Testers (parents)", nt - 1, st)
      rows$plt <- .bp_rows("    Lines vs testers", 1, ss_p - sl - st)
    }
  }
  if (np >= 1) rows$pvc <- .bp_rows("  Parents vs crosses", 1, sq$ss[["Entries"]] - ss_cross -
                                      if (np >= 2) r * sum((pm - mean(pm))^2) else 0)
  rows$c <- .bp_rows("  Crosses", l * t - 1, ss_cross)
  rows$l <- .bp_rows("    Lines", l - 1, ss_l)
  rows$t <- .bp_rows("    Testers", t - 1, ss_t)
  rows$lt <- .bp_rows("    Lines x Testers", (l - 1) * (t - 1), ss_lt)
  rows$err <- .bp_rows("Error", dfe, mse * dfe)
  rows$tot <- .bp_rows("Total", length(y) - 1, sum((y - mean(y))^2))
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  tab <- .bp_against(tab, setdiff(tab$source, c("Error", "Total")), "Error")
  lxt <- .bp_rows(c("Lines", "Testers", "Lines x Testers", "Error"), c(l - 1, t - 1, (l - 1) * (t - 1), dfe),
                  c(ss_l, ss_t, ss_lt, mse * dfe))
  lxt <- .bp_against(lxt, c("Lines", "Testers"), "Lines x Testers")
  lxt <- .bp_against(lxt, "Lines x Testers", "Error")
  # effects
  se <- c(gca_line = sqrt(mse / (r * t)), gca_tester = sqrt(mse / (r * l)), sca = sqrt(mse / r),
          diff_gca_line = sqrt(2 * mse / (r * t)), diff_gca_tester = sqrt(2 * mse / (r * l)),
          diff_sca = sqrt(2 * mse / r))
  gl <- lm_ - gm; gt <- tm - gm
  S <- M - outer(lm_, tm, "+") + gm
  sl <- .bp_signif(gl, se[["gca_line"]], dfe); stt <- .bp_signif(gt, se[["gca_tester"]], dfe)
  gca_l <- data.frame(line = L, mean = unname(lm_), gca = unname(gl), se = se[["gca_line"]], t = unname(sl$t),
                      p = unname(sl$p), stringsAsFactors = FALSE)
  gca_t <- data.frame(tester = Tt, mean = unname(tm), gca = unname(gt), se = se[["gca_tester"]], t = unname(stt$t),
                      p = unname(stt$p), stringsAsFactors = FALSE)
  ss_ <- .bp_signif(as.vector(S), se[["sca"]], dfe)
  sca <- data.frame(line = rep(L, t), tester = rep(Tt, each = l), mean = as.vector(M), sca = as.vector(S),
                    se = se[["sca"]], t = ss_$t, p = ss_$p, stringsAsFactors = FALSE)
  sca <- sca[order(match(sca$line, L), match(sca$tester, Tt)), ]; rownames(sca) <- NULL
  # genetic components (Kempthorne 1957)
  Ml <- ss_l / (l - 1); Mt <- ss_t / (t - 1); Mlt <- ss_lt / ((l - 1) * (t - 1))
  cov_l <- (Ml - Mlt) / (r * t); cov_t <- (Mt - Mlt) / (r * l)
  cov_hs <- (((l - 1) * Ml + (t - 1) * Mt) / (l + t - 2) - Mlt) * (l + t - 2) / (r * (2 * l * t - l - t))
  cov_fs <- ((Ml - mse) + (Mt - mse) + (Mlt - mse) + 6 * r * cov_hs - r * (l + t) * cov_hs) / (3 * r)
  s2_sca <- (Mlt - mse) / r
  comp <- data.frame(component = c("Cov HS (lines)", "Cov HS (testers)", "Cov HS (average)", "Cov FS",
                                   "sigma2 GCA (= Cov HS)", "sigma2 SCA"),
                     estimate = c(cov_l, cov_t, cov_hs, cov_fs, cov_hs, s2_sca), stringsAsFactors = FALSE)
  gen <- data.frame(F = c(0, 1), sigma2_A = 4 * cov_hs / (1 + c(0, 1)), sigma2_D = 4 * s2_sca / (1 + c(0, 1))^2)
  gen$dominance <- ifelse(gen$sigma2_A > 0 & gen$sigma2_D >= 0, sqrt(2 * gen$sigma2_D / gen$sigma2_A), NA_real_)
  contrib <- c(lines = 100 * ss_l / ss_cross, testers = 100 * ss_t / ss_cross, lxt = 100 * ss_lt / ss_cross)
  structure(list(anova = tab, anova_lxt = lxt, gca_lines = gca_l, gca_testers = gca_t, sca = sca, se = se,
                 contribution = contrib, components = comp, genetic = gen, means = M,
                 parents = if (np) data.frame(parent = sub("^[^\r]*\r", "", names(pm)), type = pk, mean = unname(pm),
                                              stringsAsFactors = FALSE) else NULL,
                 r = r, mse = mse, df_error = dfe, response = response, n = length(y),
                 checks_dropped = sum(checks)),
            class = "ad_line_tester")
}

#' @export
print.ad_line_tester <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  l <- nrow(x$means); t <- ncol(x$means)
  cat("\nLine x tester analysis (Kempthorne 1957)\n")
  cat("Response: ", x$response, "   ", l, " lines x ", t, " testers, ",
      if (is.null(x$parents)) 0 else nrow(x$parents), " parents, r = ", x$r,
      if (x$checks_dropped) paste0("; ", x$checks_dropped, " row(s) without line or tester not used"), "\n\n", sep = "")
  .bp_print_anova(x$anova, digits)
  cat("\nLine x tester analysis (random model: lines and testers against L x T):\n")
  .bp_print_anova(x$anova_lxt, digits)
  cat("\nGCA effects of lines:\n"); .bp_print_df(x$gca_lines, digits)
  cat("\nGCA effects of testers:\n"); .bp_print_df(x$gca_testers, digits)
  cat("\nSCA effects:\n")
  S <- x$means * 0
  S[cbind(match(x$sca$line, rownames(S)), match(x$sca$tester, colnames(S)))] <- x$sca$sca
  print(round(S, max(2, digits - 2)))
  cat("\nStandard errors:\n")
  se <- x$se
  names(se) <- c("GCA line", "GCA tester", "SCA", "GCA line difference", "GCA tester difference", "SCA difference")
  print(signif(se, digits))
  cat("\nGenetic components:\n"); .bp_print_df(x$components, digits)
  cat("\n"); .bp_print_df(x$genetic, digits)
  cat("\nContribution to the variation among crosses (%):\n")
  print(round(x$contribution, 2))
  invisible(x)
}

# ---------- Griffing's diallels ----------

# Parameters: mu, g (p), s (pairs), m (p, maternal), n (reciprocal pairs, non-maternal).
.dl_layout <- function(p, method) {
  selfs <- method %in% 1:2; recip <- method %in% c(1, 3)
  pairs <- NULL
  for (i in seq_len(p)) for (j in i:p) if (i < j || selfs) pairs <- rbind(pairs, c(i, j))
  pix <- matrix(NA_integer_, p, p)
  for (k in seq_len(nrow(pairs))) pix[pairs[k, 1], pairs[k, 2]] <- pix[pairs[k, 2], pairs[k, 1]] <- k
  rpairs <- NULL
  if (recip) for (i in seq_len(p - 1)) for (j in (i + 1):p) rpairs <- rbind(rpairs, c(i, j))
  rix <- matrix(NA_integer_, p, p)
  if (recip) for (k in seq_len(nrow(rpairs))) rix[rpairs[k, 1], rpairs[k, 2]] <- rix[rpairs[k, 2], rpairs[k, 1]] <- k
  off <- c(mu = 1, g = 1, s = 1 + p, m = 1 + p + nrow(pairs), n = 1 + 2 * p + nrow(pairs))
  npar <- 1 + p + nrow(pairs) + if (recip) p + nrow(rpairs) else 0
  list(p = p, method = method, selfs = selfs, recip = recip, pairs = pairs, pix = pix, rpairs = rpairs,
       rix = rix, off = off, npar = npar)
}

.dl_design <- function(Lt, i, j) {
  X <- matrix(0, length(i), Lt$npar)
  for (k in seq_along(i)) {
    a <- i[k]; b <- j[k]
    X[k, 1] <- 1
    X[k, Lt$off[["g"]] + a] <- X[k, Lt$off[["g"]] + a] + 1
    X[k, Lt$off[["g"]] + b] <- X[k, Lt$off[["g"]] + b] + 1
    X[k, Lt$off[["s"]] + Lt$pix[a, b]] <- 1
    if (Lt$recip && a != b) {
      X[k, Lt$off[["m"]] + a] <- 1; X[k, Lt$off[["m"]] + b] <- -1
      X[k, Lt$off[["n"]] + Lt$rix[a, b]] <- if (a < b) 1 else -1
    }
  }
  X
}

.dl_restrictions <- function(Lt) {
  p <- Lt$p; C <- NULL
  z <- function() numeric(Lt$npar)
  r <- z(); r[Lt$off[["g"]] + seq_len(p)] <- 1; C <- rbind(C, r)
  for (i in seq_len(p)) {
    r <- z()
    for (j in seq_len(p)) {
      if (i == j) { if (Lt$selfs) r[Lt$off[["s"]] + Lt$pix[i, i]] <- r[Lt$off[["s"]] + Lt$pix[i, i]] + if (Lt$method == 2) 2 else 1 }
      else r[Lt$off[["s"]] + Lt$pix[i, j]] <- r[Lt$off[["s"]] + Lt$pix[i, j]] + 1
    }
    C <- rbind(C, r)
  }
  if (Lt$recip) {
    r <- z(); r[Lt$off[["m"]] + seq_len(p)] <- 1; C <- rbind(C, r)
    for (i in seq_len(p)) {
      r <- z()
      for (j in seq_len(p)[-i]) r[Lt$off[["n"]] + Lt$rix[i, j]] <- if (i < j) 1 else -1
      C <- rbind(C, r)
    }
  }
  C
}

# Least squares subject to C theta = 0, on a basis of the null space of C.
.bp_cls <- function(X, y, C) {
  qc <- qr(t(C))
  N <- qr.Q(qc, complete = TRUE)[, -seq_len(qc$rank), drop = FALSE]
  Z <- X %*% N
  Zi <- solve(crossprod(Z))
  b <- Zi %*% crossprod(Z, y)
  list(theta = drop(N %*% b), V = N %*% Zi %*% t(N))
}

#' Diallel analysis by Griffing's methods
#'
#' Griffing's (1956) analysis of a diallel cross by any of the four methods:
#' 1, parents, F1 and reciprocals; 2, parents and one F1 per pair; 3, F1 and
#' reciprocals; 4, one F1 per pair. The analysis of variance of general
#' (GCA) and specific (SCA) combining ability, and of reciprocal effects
#' (split into maternal and non-maternal, Cockerham 1963) in methods 1 and
#' 3; the effects with their exact standard errors (Model I), the standard
#' errors of differences, and the variance components of Model II with the
#' additive and dominance variances.
#'
#' The combining-ability analysis is on the scale of entry means: the error
#' is the error mean square of the trial divided by the number of
#' replicates. Rows that do not belong to the method are left out (parents
#' in methods 3 and 4; reciprocals in methods 2 and 4, where the cross with
#' the first parent as female is kept), so one data set serves all methods.
#' The effects are the least-squares estimates under Griffing's
#' restrictions (sum of GCA effects zero; SCA effects summing to zero for
#' every parent, the selfs counted twice in method 2), which reproduce his
#' closed formulas.
#'
#' @param data A data frame in long format.
#' @param response Name of the numeric response column.
#' @param parent1,parent2 Names of the female and male parent columns
#'   (equal for the parents themselves).
#' @param block Name of the block column (randomised complete blocks); when
#'   `NULL` the entries are taken as completely randomised.
#' @param method Griffing's method, 1 to 4.
#' @param inbreeding Inbreeding coefficient of the parents, used to turn the
#'   Model II components into additive and dominance variances (1 for inbred
#'   lines).
#' @return An object of class `"ad_diallel"`.
#' @references Griffing, B. (1956) Australian Journal of Biological Sciences
#'   9, 463-493. Cockerham, C. C. (1963) In Statistical Genetics and Plant
#'   Breeding, NAS-NRC 982, 53-94.
#' @examples
#' dl <- diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block", method = 1)
#' dl
#' diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block", method = 4)$gca
#' @export
diallel <- function(data, response, parent1, parent2, block = NULL, method = 1, inbreeding = 1) {
  if (!is.numeric(method) || length(method) != 1 || !method %in% 1:4)
    stop("'method' must be 1, 2, 3 or 4 (Griffing 1956)", call. = FALSE)
  .bp_check_inbreeding(inbreeding)
  d <- .bp_input(data, response, c(parent1, parent2, block))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  a <- d[[parent1]]; b <- d[[parent2]]
  lv <- if (is.factor(a) || is.factor(b)) unique(c(levels(as.factor(a)), levels(as.factor(b))))
        else .levels_of(c(a, b))
  lv <- lv[lv %in% c(as.character(a), as.character(b))]
  i <- match(as.character(a), lv); j <- match(as.character(b), lv)
  p <- length(lv)
  Lt <- .dl_layout(p, method)
  keep <- rep(TRUE, length(i))
  if (!Lt$selfs) keep[i == j] <- FALSE
  if (!Lt$recip) {
    key <- paste(pmin(i, j), pmax(i, j))
    direct <- tapply(i <= j, key, any)
    keep[i > j & direct[key]] <- FALSE
    sw <- keep & i > j
    tmp <- i[sw]; i[sw] <- j[sw]; j[sw] <- tmp
  }
  n_left <- sum(!keep)
  d <- d[keep, , drop = FALSE]; i <- i[keep]; j <- j[keep]
  if (p < 3) stop("at least three parents are needed", call. = FALSE)
  if (method %in% 3:4 && p < 4) stop("methods 3 and 4 need at least four parents", call. = FALSE)
  # every cell of the method
  cells <- if (Lt$recip) {
    cc <- expand.grid(i = seq_len(p), j = seq_len(p)); cc[Lt$selfs | cc$i != cc$j, ]
  } else data.frame(i = Lt$pairs[, 1], j = Lt$pairs[, 2])
  ck <- paste(cells$i, cells$j); dk <- paste(i, j)
  miss <- !ck %in% dk
  if (any(miss)) stop("crosses missing for method ", method, ": ",
                      paste(utils::head(paste(lv[cells$i[miss]], lv[cells$j[miss]], sep = " x "), 5), collapse = ", "),
                      if (sum(miss) > 5) ", ..." else "", call. = FALSE)
  Ent <- factor(dk, levels = ck)
  y <- d[[response]]
  if (!is.null(block)) {
    B <- factor(as.character(d[[block]]))
    if (nlevels(B) < 2) stop("at least two blocks are needed", call. = FALSE)
    if (any(table(B, Ent) != 1)) stop("every entry must appear exactly once in every block", call. = FALSE)
    r <- nlevels(B)
    sq <- .bp_seq(y, list(Blocks = B, Entries = Ent))
  } else {
    cnt <- table(Ent)
    if (length(unique(cnt)) != 1 || cnt[1] < 2)
      stop("every entry needs the same number of replicates (at least two)", call. = FALSE)
    r <- unname(cnt[1])
    sq <- .bp_seq(y, list(Entries = Ent))
  }
  dfe <- sq$df_res; mse <- sq$rss / dfe
  if (dfe < 1) stop("no degrees of freedom left for the error", call. = FALSE)
  rows <- list()
  if (!is.null(block)) rows$b <- .bp_rows("Blocks", sq$df[["Blocks"]], sq$ss[["Blocks"]])
  rows$e <- .bp_rows("Entries", sq$df[["Entries"]], sq$ss[["Entries"]])
  rows$err <- .bp_rows("Error", dfe, mse * dfe)
  rows$tot <- .bp_rows("Total", length(y) - 1, sum((y - mean(y))^2))
  tab_e <- do.call(rbind, rows); rownames(tab_e) <- NULL
  tab_e <- .bp_against(tab_e, c("Blocks", "Entries"), "Error")
  # combining ability on the scale of entry means
  ym <- tapply(y, Ent, mean)
  ci <- cells$i; cj <- cells$j
  Gi <- matrix(0, length(ci), p); for (k in seq_along(ci)) { Gi[k, ci[k]] <- Gi[k, ci[k]] + 1; Gi[k, cj[k]] <- Gi[k, cj[k]] + 1 }
  Si <- .bp_ind(factor(Lt$pix[cbind(ci, cj)], levels = seq_len(nrow(Lt$pairs))))
  terms <- list(GCA = Gi, SCA = Si)
  if (Lt$recip) {
    Mi <- matrix(0, length(ci), p); Ni <- matrix(0, length(ci), nrow(Lt$rpairs))
    for (k in seq_along(ci)) if (ci[k] != cj[k]) {
      Mi[k, ci[k]] <- 1; Mi[k, cj[k]] <- -1
      Ni[k, Lt$rix[ci[k], cj[k]]] <- if (ci[k] < cj[k]) 1 else -1
    }
    terms$MAT <- Mi; terms$NONMAT <- Ni
  }
  sc <- .bp_seq(as.vector(ym), terms)
  mse_m <- mse / r
  ca <- .bp_rows(c("GCA", "SCA"), sc$df[c("GCA", "SCA")], sc$ss[c("GCA", "SCA")])
  if (Lt$recip)
    ca <- rbind(ca, .bp_rows(c("Reciprocal", "  Maternal", "  Non-maternal"),
                             c(sc$df[["MAT"]] + sc$df[["NONMAT"]], sc$df[["MAT"]], sc$df[["NONMAT"]]),
                             c(sc$ss[["MAT"]] + sc$ss[["NONMAT"]], sc$ss[["MAT"]], sc$ss[["NONMAT"]])))
  ca <- rbind(ca, .bp_rows("Error", dfe, mse_m * dfe))
  ca <- .bp_against(ca, setdiff(ca$source, "Error"), "Error")
  ms <- stats::setNames(ca$ms, trimws(ca$source))
  F2 <- ms[["GCA"]] / ms[["SCA"]]
  gca_ii <- list(F = F2, df1 = ca$df[1], df2 = ca$df[2],
                 p = stats::pf(F2, ca$df[1], ca$df[2], lower.tail = FALSE), exact = method != 1)
  # effects (Model I)
  fit <- .bp_cls(.dl_design(Lt, ci, cj), as.vector(ym), .dl_restrictions(Lt))
  th <- fit$theta; V <- fit$V * mse_m
  eff <- function(idx, ...) {
    est <- th[idx]; se <- sqrt(pmax(0, diag(V)[idx])); s <- .bp_signif(est, se, dfe)
    data.frame(..., estimate = est, se = se, t = s$t, p = s$p, stringsAsFactors = FALSE, row.names = NULL)
  }
  gca <- eff(Lt$off[["g"]] + seq_len(p), parent = lv)
  sca <- eff(Lt$off[["s"]] + seq_len(nrow(Lt$pairs)), parent1 = lv[Lt$pairs[, 1]], parent2 = lv[Lt$pairs[, 2]])
  vdiff <- function(a, b) V[cbind(a, a)] + V[cbind(b, b)] - 2 * V[cbind(a, b)]
  pr <- t(utils::combn(p, 2))
  sed <- data.frame(difference = "g_i - g_j", se = sqrt(mean(vdiff(Lt$off[["g"]] + pr[, 1], Lt$off[["g"]] + pr[, 2]))),
                    stringsAsFactors = FALSE)
  sp <- Lt$pairs; np <- nrow(sp)
  if (np > 1) {
    pp <- t(utils::combn(np, 2))
    self1 <- sp[pp[, 1], 1] == sp[pp[, 1], 2]; self2 <- sp[pp[, 2], 1] == sp[pp[, 2], 2]
    share <- apply(pp, 1, function(z) length(intersect(sp[z[1], ], sp[z[2], ])) > 0)
    vs <- vdiff(Lt$off[["s"]] + pp[, 1], Lt$off[["s"]] + pp[, 2])
    add <- function(lab, sel) if (any(sel)) sed[nrow(sed) + 1, ] <<- list(lab, sqrt(mean(vs[sel])))
    add("s_ii - s_jj", self1 & self2)
    add("s_ij - s_ik", !self1 & !self2 & share)
    add("s_ij - s_kl", !self1 & !self2 & !share)
  }
  out <- list(method = method, parents = lv, p = p, anova_entries = tab_e, anova = ca, gca = gca, sca = sca)
  if (Lt$recip) {
    nr <- nrow(Lt$rpairs)
    Lr <- matrix(0, nr, Lt$npar)
    for (k in seq_len(nr)) {
      a <- Lt$rpairs[k, 1]; bb <- Lt$rpairs[k, 2]
      Lr[k, Lt$off[["m"]] + a] <- 1; Lr[k, Lt$off[["m"]] + bb] <- -1; Lr[k, Lt$off[["n"]] + k] <- 1
    }
    est <- drop(Lr %*% th); vr <- diag(Lr %*% V %*% t(Lr)); s <- .bp_signif(est, sqrt(vr), dfe)
    out$reciprocal <- data.frame(parent1 = lv[Lt$rpairs[, 1]], parent2 = lv[Lt$rpairs[, 2]], estimate = est,
                                 se = sqrt(vr), t = s$t, p = s$p, stringsAsFactors = FALSE)
    out$maternal <- eff(Lt$off[["m"]] + seq_len(p), parent = lv)
    if (nr > 1) {
      rr <- t(utils::combn(nr, 2))
      vd <- vapply(seq_len(nrow(rr)), function(k) { w <- Lr[rr[k, 1], ] - Lr[rr[k, 2], ]; sum(w * (V %*% w)) }, 0)
      sed[nrow(sed) + 1, ] <- list("r_ij - r_kl", sqrt(mean(vd)))
    }
  }
  out$se_differences <- sed
  # Model I: quadratic forms of the fixed effects
  cg <- switch(method, 2 * p, p + 2, 2 * (p - 2), p - 2)
  cs <- if (method %in% c(1, 3)) 2 else 1
  m1 <- data.frame(quantity = c("sum(g^2)/(p - 1)", "sum(s^2)/df"),
                   estimate = c((ms[["GCA"]] - mse_m) / cg, (ms[["SCA"]] - mse_m) / cs), stringsAsFactors = FALSE)
  if (Lt$recip) m1[3, ] <- list("sum(r^2)/df", (ms[["Reciprocal"]] - mse_m) / 2)
  out$model1 <- m1
  # Model II: variance components (Griffing's expected mean squares)
  dfs <- stats::setNames(ca$df, trimws(ca$source))
  if (method == 1) {
    den <- p^2 - p + 1
    s2s <- p^2 * (ms[["SCA"]] - mse_m) / (2 * den)
    se_s <- p^2 / (2 * den) * .bp_lin_se(c(1, -1), c(ms[["SCA"]], mse_m), c(dfs[["SCA"]], dfe))
    s2g <- (ms[["GCA"]] - (mse_m + p * (p - 1) * ms[["SCA"]]) / den) / (2 * p)
    se_g <- .bp_lin_se(c(1 / (2 * p), -1 / (2 * p * den), -p * (p - 1) / (2 * p * den)),
                       c(ms[["GCA"]], mse_m, ms[["SCA"]]), c(dfs[["GCA"]], dfe, dfs[["SCA"]]))
  } else {
    css <- if (method == 3) 2 else 1
    s2s <- (ms[["SCA"]] - mse_m) / css
    se_s <- .bp_lin_se(c(1, -1) / css, c(ms[["SCA"]], mse_m), c(dfs[["SCA"]], dfe))
    s2g <- (ms[["GCA"]] - ms[["SCA"]]) / cg
    se_g <- .bp_lin_se(c(1, -1) / cg, c(ms[["GCA"]], ms[["SCA"]]), c(dfs[["GCA"]], dfs[["SCA"]]))
  }
  m2 <- data.frame(component = c("sigma2 GCA", "sigma2 SCA"), estimate = c(s2g, s2s), se = c(se_g, se_s),
                   stringsAsFactors = FALSE)
  if (Lt$recip) {
    m2[3, ] <- list("sigma2 reciprocal", (ms[["Reciprocal"]] - mse_m) / 2,
                    .bp_lin_se(c(0.5, -0.5), c(ms[["Reciprocal"]], mse_m), c(dfs[["Reciprocal"]], dfe)))
    m2[4, ] <- list("sigma2 maternal", (ms[["Maternal"]] - ms[["Non-maternal"]]) / (2 * p),
                    .bp_lin_se(c(1, -1) / (2 * p), c(ms[["Maternal"]], ms[["Non-maternal"]]),
                               c(dfs[["Maternal"]], dfs[["Non-maternal"]])))
    m2[5, ] <- list("sigma2 non-maternal", (ms[["Non-maternal"]] - mse_m) / 2,
                    .bp_lin_se(c(0.5, -0.5), c(ms[["Non-maternal"]], mse_m), c(dfs[["Non-maternal"]], dfe)))
  }
  m2[nrow(m2) + 1, ] <- list("sigma2 error (entry mean)", mse_m, .bp_lin_se(1, mse_m, dfe))
  out$model2 <- m2
  Fi <- inbreeding
  s2A <- 4 * s2g / (1 + Fi); s2D <- 4 * s2s / (1 + Fi)^2
  out$genetic <- data.frame(parameter = c("sigma2_A", "sigma2_D", "dominance ratio sqrt(2 sigma2_D / sigma2_A)",
                                          "Baker ratio 2 sigma2_GCA / (2 sigma2_GCA + sigma2_SCA)"),
                            estimate = c(s2A, s2D, if (s2A > 0 && s2D >= 0) sqrt(2 * s2D / s2A) else NA_real_,
                                         if (s2g > 0 && s2s >= 0) 2 * s2g / (2 * s2g + s2s) else NA_real_),
                            stringsAsFactors = FALSE)
  out$gca_model2 <- gca_ii
  out$means <- data.frame(parent1 = lv[ci], parent2 = lv[cj], mean = as.vector(ym), stringsAsFactors = FALSE)
  out$inbreeding <- inbreeding; out$r <- r; out$mse <- mse; out$df_error <- dfe
  out$response <- response; out$n <- length(y); out$rows_left_out <- n_left
  structure(out, class = "ad_diallel")
}

#' @export
print.ad_diallel <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  lab <- c("parents, F1 and reciprocals", "parents and F1", "F1 and reciprocals", "F1 only")[x$method]
  cat("\nDiallel analysis, Griffing (1956) method ", x$method, " (", lab, ")\n", sep = "")
  cat("Response: ", x$response, "   ", x$p, " parents, r = ", x$r,
      if (x$rows_left_out) paste0("; ", x$rows_left_out, " row(s) outside this method left out"), "\n\n", sep = "")
  cat("Analysis of variance of the entries (plot scale):\n")
  .bp_print_anova(x$anova_entries, digits)
  cat("\nCombining ability (scale of entry means; error = MSE / r):\n")
  .bp_print_anova(x$anova, digits)
  g <- x$gca_model2
  cat(sprintf("Model II: GCA against SCA, F = %s on %d and %d df, p = %s%s\n", format(g$F, digits = digits),
              as.integer(g$df1), as.integer(g$df2), format.pval(g$p, digits = 3),
              if (g$exact) "" else " (approximate in method 1)"))
  cat("\nGCA effects:\n"); .bp_print_df(x$gca, digits)
  cat("\nSCA effects:\n"); .bp_print_df(x$sca, digits)
  if (!is.null(x$reciprocal)) {
    cat("\nReciprocal effects r_ij = (x_ij - x_ji) / 2:\n"); .bp_print_df(x$reciprocal, digits)
    cat("\nMaternal effects:\n"); .bp_print_df(x$maternal, digits)
  }
  cat("\nStandard errors of differences:\n"); .bp_print_df(x$se_differences, digits)
  cat("\nModel I (fixed effects):\n"); .bp_print_df(x$model1, digits)
  cat("\nModel II (random effects), entry-mean scale:\n"); .bp_print_df(x$model2, digits)
  cat("\nGenetic parameters (inbreeding of the parents F = ", x$inbreeding, "):\n", sep = "")
  .bp_print_df(x$genetic, digits)
  invisible(x)
}

# ---------- North Carolina designs ----------

#' North Carolina designs I, II and III
#'
#' Analysis of the mating designs of Comstock and Robinson (1948, 1952):
#' \describe{
#'   \item{Design I}{each male crossed with its own group of females
#'     (females nested within males).}
#'   \item{Design II}{every male crossed with every female of its set.}
#'   \item{Design III}{F2 (or other) plants used as males, each crossed with
#'     the two inbred parents of the population (the testers, given in the
#'     `female` column).}
#' }
#' Males (and females) are nested within sets, and blocks within sets. When
#' a plot holds several plants, the plot error (blocks x crosses within
#' sets) is separated from the variation within plots. The table gives the
#' expected mean squares; every term is tested against the mean square that
#' differs from it only by its own component, and the components are their
#' moment estimators with standard errors from Var(MS) = 2 MS^2 / df.
#'
#' Genetic interpretation without epistasis: design I,
#' `sigma2_A = 4 Vm / (1 + F)` and `sigma2_D = 4 (Vf - Vm) / (1 + F)^2`;
#' design II, `sigma2_A = 2 (Vm + Vf) / (1 + F)` and
#' `sigma2_D = 4 Vmf / (1 + F)^2`; design III (F2 reference population,
#' gene frequency 1/2), `sigma2_A = 4 Vm = 2 (MSm - MSe) / r` and
#' `sigma2_D = Vmt = (MSmt - MSe) / r`, because the sum and the difference
#' of the two backcross families of a plant have variances sigma2_A and
#' 2 sigma2_D. The average degree of dominance is `sqrt(2 sigma2_D / sigma2_A)`.
#'
#' @param data A data frame in long format.
#' @param response Name of the numeric response column.
#' @param male,female Names of the male and female columns (design III: the
#'   female column holds the two testers).
#' @param block Name of the block (replicate) column, nested within sets.
#' @param set Optional name of the set column.
#' @param design 1, 2 or 3.
#' @param inbreeding Inbreeding coefficient of the parents (designs I and II).
#' @return An object of class `"ad_nc_design"`: the analysis of variance with
#'   expected mean squares, the variance components and the genetic
#'   parameters.
#' @references Comstock, R. E. and Robinson, H. F. (1948) Biometrics 4,
#'   254-266. Comstock, R. E. and Robinson, H. F. (1952) In Heterosis, Iowa
#'   State College Press, 494-516.
#' @examples
#' nc_design(nc_designs$nc2, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 2)
#' nc_design(nc_designs$nc3, "Yield_t_ha", "Male", "Tester", block = "Rep", set = "Set", design = 3)
#' @export
nc_design <- function(data, response, male, female, block, set = NULL, design = 1, inbreeding = 0) {
  if (!is.numeric(design) || length(design) != 1 || !design %in% 1:3)
    stop("'design' must be 1, 2 or 3 (North Carolina design I, II or III)", call. = FALSE)
  .bp_check_inbreeding(inbreeding)
  d <- .bp_input(data, response, c(male, female, block, set))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  y <- d[[response]]
  S <- if (is.null(set)) factor(rep("1", nrow(d))) else factor(as.character(d[[set]]))
  s <- nlevels(S)
  R <- .bp_fac(S, d[[block]])
  Ml <- .bp_fac(S, d[[male]])
  Tt <- NULL
  if (design == 1) { Fe <- .bp_fac(S, d[[male]], d[[female]]); Cr <- Fe }
  if (design == 2) { Fe <- .bp_fac(S, d[[female]]); Cr <- .bp_fac(S, d[[male]], d[[female]]) }
  if (design == 3) {
    Tt <- factor(as.character(d[[female]]))
    if (nlevels(Tt) != 2) stop("design III needs exactly two testers (the inbred parents) in '", female, "'", call. = FALSE)
    Cr <- .bp_fac(S, d[[male]], d[[female]])
  }
  # balance
  set_of <- function(f) S[match(levels(f), f)]
  per <- function(f, g) { u <- unique(data.frame(f = f, g = g)); table(u$g) }
  rs <- per(R, S); ms_ <- per(Ml, S)
  if (length(unique(rs)) != 1 || rs[1] < 2) stop("every set needs the same number of blocks (at least two)", call. = FALSE)
  if (length(unique(ms_)) != 1 || ms_[1] < 2) stop("every set needs the same number of males (at least two)", call. = FALSE)
  r <- unname(rs[1]); m <- unname(ms_[1])
  if (design == 1) {
    fm <- per(Fe, Ml)
    if (length(unique(fm)) != 1 || fm[1] < 2) stop("every male needs the same number of females (at least two)", call. = FALSE)
    f <- unname(fm[1])
  } else if (design == 2) {
    fs <- per(Fe, S)
    if (length(unique(fs)) != 1 || fs[1] < 2) stop("every set needs the same number of females (at least two)", call. = FALSE)
    f <- unname(fs[1])
    if (any(per(Cr, S) != m * f)) stop("design II needs every male crossed with every female of its set", call. = FALSE)
  } else {
    f <- 2
    if (any(per(Cr, Ml) != 2)) stop("design III needs every male crossed with both testers", call. = FALSE)
  }
  tc <- table(Cr, R)
  P <- .bp_fac(R, Cr)
  npl <- table(P)
  if (length(unique(npl)) != 1) stop("every plot must hold the same number of observations", call. = FALSE)
  k <- unname(npl[1])
  if (any(tc != outer(as.character(set_of(Cr)), as.character(set_of(R)), "==") * k))
    stop("every cross must appear once in every block of its set", call. = FALSE)
  # analysis of variance
  terms <- list()
  if (s > 1) terms$Sets <- S
  terms$`Blocks(sets)` <- R
  if (design == 1) { terms$`Males(sets)` <- Ml; terms$`Females(males)` <- Fe }
  if (design == 2) { terms$`Males(sets)` <- Ml; terms$`Females(sets)` <- Fe; terms$`Males x Females(sets)` <- Cr }
  if (design == 3) {
    terms$Testers <- Tt
    if (s > 1) terms$`Testers x Sets` <- .bp_fac(Tt, S)
    terms$`Males(sets)` <- Ml; terms$`Males x Testers(sets)` <- Cr
  }
  if (k > 1) terms$`Plot error` <- P
  sq <- .bp_seq(y, terms)
  res_name <- if (k > 1) "Within plots" else "Error"
  tab <- rbind(.bp_rows(names(terms), unname(sq$df), unname(sq$ss)), .bp_rows(res_name, sq$df_res, sq$rss))
  tab <- rbind(tab, .bp_rows("Total", length(y) - 1, sum((y - mean(y))^2)))
  err <- if (k > 1) "Plot error" else "Error"
  ms <- stats::setNames(tab$ms, tab$source); dfs <- stats::setNames(tab$df, tab$source)
  base <- if (k > 1) sprintf("Vw + %d Ve", k) else "Ve"
  nr <- k * r
  ems <- stats::setNames(rep("", nrow(tab)), tab$source)
  ems[err] <- base
  if (k > 1) ems["Within plots"] <- "Vw"
  if (design == 1) {
    ems["Females(males)"] <- sprintf("%s + %d Vf", base, nr)
    ems["Males(sets)"] <- sprintf("%s + %d Vf + %d Vm", base, nr, nr * f)
    tab <- .bp_against(tab, "Males(sets)", "Females(males)")
    tab <- .bp_against(tab, c("Females(males)", "Blocks(sets)"), err)
    Vm <- (ms[["Males(sets)"]] - ms[["Females(males)"]]) / (nr * f)
    Vf <- (ms[["Females(males)"]] - ms[[err]]) / nr
    cm <- c(1, -1) / (nr * f); cf <- c(1, -1) / nr
    comp <- data.frame(component = c("Vm (males)", "Vf (females within males)"), estimate = c(Vm, Vf),
                       se = c(.bp_lin_se(cm, ms[c("Males(sets)", "Females(males)")], dfs[c("Males(sets)", "Females(males)")]),
                              .bp_lin_se(cf, ms[c("Females(males)", err)], dfs[c("Females(males)", err)])),
                       stringsAsFactors = FALSE)
    s2A <- 4 * Vm / (1 + inbreeding); s2D <- 4 * (Vf - Vm) / (1 + inbreeding)^2
    # Vm + Vf and Vf - Vm as combinations of the mean squares of males, females and error
    three <- c("Males(sets)", "Females(males)", err)
    se_FS <- .bp_lin_se(c(1 / (nr * f), 1 / nr - 1 / (nr * f), -1 / nr), ms[three], dfs[three])
    se_D <- 4 / (1 + inbreeding)^2 * .bp_lin_se(c(-1 / (nr * f), 1 / nr + 1 / (nr * f), -1 / nr), ms[three], dfs[three])
    gen <- data.frame(parameter = c("Cov HS = Vm", "Cov FS = Vm + Vf", "sigma2_A", "sigma2_D"),
                      estimate = c(Vm, Vm + Vf, s2A, s2D),
                      se = c(comp$se[1], se_FS, 4 / (1 + inbreeding) * comp$se[1], se_D),
                      stringsAsFactors = FALSE)
  } else if (design == 2) {
    ems["Males x Females(sets)"] <- sprintf("%s + %d Vmf", base, nr)
    ems["Males(sets)"] <- sprintf("%s + %d Vmf + %d Vm", base, nr, nr * f)
    ems["Females(sets)"] <- sprintf("%s + %d Vmf + %d Vf", base, nr, nr * m)
    tab <- .bp_against(tab, c("Males(sets)", "Females(sets)"), "Males x Females(sets)")
    tab <- .bp_against(tab, c("Males x Females(sets)", "Blocks(sets)"), err)
    mf <- "Males x Females(sets)"
    Vm <- (ms[["Males(sets)"]] - ms[[mf]]) / (nr * f)
    Vf <- (ms[["Females(sets)"]] - ms[[mf]]) / (nr * m)
    Vmf <- (ms[[mf]] - ms[[err]]) / nr
    se3 <- c(.bp_lin_se(c(1, -1) / (nr * f), ms[c("Males(sets)", mf)], dfs[c("Males(sets)", mf)]),
             .bp_lin_se(c(1, -1) / (nr * m), ms[c("Females(sets)", mf)], dfs[c("Females(sets)", mf)]),
             .bp_lin_se(c(1, -1) / nr, ms[c(mf, err)], dfs[c(mf, err)]))
    comp <- data.frame(component = c("Vm (males)", "Vf (females)", "Vmf (males x females)"),
                       estimate = c(Vm, Vf, Vmf), se = se3, stringsAsFactors = FALSE)
    a <- 4 / (1 + inbreeding)
    # (Vm + Vf) / 2 as a combination of the mean squares
    cA <- c(1 / (2 * nr * f), 1 / (2 * nr * m), -1 / (2 * nr * f) - 1 / (2 * nr * m))
    gen <- data.frame(parameter = c("sigma2_A (males)", "sigma2_A (females)", "sigma2_A (pooled)", "sigma2_D"),
                      estimate = c(a * Vm, a * Vf, a * (Vm + Vf) / 2, 4 * Vmf / (1 + inbreeding)^2),
                      se = c(a * se3[1], a * se3[2],
                             a * .bp_lin_se(cA, ms[c("Males(sets)", "Females(sets)", mf)], dfs[c("Males(sets)", "Females(sets)", mf)]),
                             4 / (1 + inbreeding)^2 * se3[3]),
                      stringsAsFactors = FALSE)
    s2A <- gen$estimate[3]; s2D <- gen$estimate[4]
  } else {
    mt <- "Males x Testers(sets)"
    ems[mt] <- sprintf("%s + %d Vmt", base, nr)
    ems["Males(sets)"] <- sprintf("%s + %d Vm", base, 2 * nr)
    ems["Testers"] <- sprintf("%s + %d Vmt + fixed testers", base, nr)
    tab <- .bp_against(tab, c("Testers", "Testers x Sets"), mt)
    tab <- .bp_against(tab, c("Males(sets)", mt, "Blocks(sets)"), err)
    Vm <- (ms[["Males(sets)"]] - ms[[err]]) / (2 * nr)
    Vmt <- (ms[[mt]] - ms[[err]]) / nr
    se2 <- c(.bp_lin_se(c(1, -1) / (2 * nr), ms[c("Males(sets)", err)], dfs[c("Males(sets)", err)]),
             .bp_lin_se(c(1, -1) / nr, ms[c(mt, err)], dfs[c(mt, err)]))
    comp <- data.frame(component = c("Vm (males)", "Vmt (males x testers)"), estimate = c(Vm, Vmt), se = se2,
                       stringsAsFactors = FALSE)
    s2A <- 4 * Vm; s2D <- Vmt
    gen <- data.frame(parameter = c("sigma2_A", "sigma2_D"), estimate = c(s2A, s2D), se = c(4 * se2[1], se2[2]),
                      stringsAsFactors = FALSE)
  }
  gen[nrow(gen) + 1, ] <- list("dominance ratio sqrt(2 sigma2_D / sigma2_A)",
                               if (s2A > 0 && s2D >= 0) sqrt(2 * s2D / s2A) else NA_real_, NA_real_)
  tab$ems <- unname(ems[tab$source])
  structure(list(design = design, anova = tab, components = comp, genetic = gen,
                 sizes = c(sets = s, blocks = r, males = m, females = f, per_plot = k),
                 inbreeding = inbreeding, response = response, n = length(y)),
            class = "ad_nc_design")
}

#' @export
print.ad_nc_design <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  nm <- c("I (females nested within males)", "II (males x females)", "III (backcrosses to both parents)")
  z <- x$sizes
  cat("\nNorth Carolina design ", nm[x$design], "\n", sep = "")
  cat("Response: ", x$response, "   sets = ", z[["sets"]], ", blocks per set = ", z[["blocks"]],
      ", males per set = ", z[["males"]],
      if (x$design == 1) paste0(", females per male = ", z[["females"]])
      else if (x$design == 2) paste0(", females per set = ", z[["females"]]) else ", testers = 2",
      if (z[["per_plot"]] > 1) paste0(", ", z[["per_plot"]], " plants per plot"), "\n\n", sep = "")
  .bp_print_anova(x$anova, digits)
  cat("\nExpected mean squares:\n")
  e <- x$anova[x$anova$ems != "", c("source", "ems")]
  for (i in seq_len(nrow(e))) cat("  ", formatC(e$source[i], width = -24), e$ems[i], "\n", sep = "")
  cat("\nVariance components:\n"); .bp_print_df(x$components, digits)
  cat("\nGenetic parameters", if (x$design < 3) paste0(" (inbreeding of the parents F = ", x$inbreeding, ")")
      else " (F2 reference population)", ":\n", sep = "")
  .bp_print_df(x$genetic, digits)
  invisible(x)
}

# ---------- example data ----------

#' Maize hybrids in six environments (multi-environment trial)
#'
#' Simulated yields of 12 maize hybrids in six environments, each a
#' randomised complete block design with three blocks; the interaction has a
#' multiplicative structure of rank two plus noise.
#'
#' @format A data frame with 216 rows and 4 columns:
#' \describe{
#'   \item{Environment}{Environment (E1-E6).}
#'   \item{Block}{Block within environment (1-3).}
#'   \item{Hybrid}{Hybrid (H01-H12).}
#'   \item{Yield_t_ha}{Grain yield, t/ha.}
#' }
#' @examples
#' ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
"maize_met"

#' Full diallel among six maize lines
#'
#' Simulated yields of a complete diallel (Griffing's method 1: parents, F1
#' and reciprocals) among six inbred lines, in three randomised complete
#' blocks. The other methods use subsets of the same data.
#'
#' @format A data frame with 108 rows and 4 columns:
#' \describe{
#'   \item{Block}{Block (1-3).}
#'   \item{Female}{Female parent (P1-P6).}
#'   \item{Male}{Male parent (P1-P6); equal to `Female` for the parents.}
#'   \item{Yield_t_ha}{Grain yield, t/ha.}
#' }
#' @examples
#' diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block", method = 2)
"maize_diallel"

#' Line x tester trial in maize
#'
#' Simulated yields of six lines crossed with three testers, evaluated with
#' the nine parents in three randomised complete blocks. Parents have an
#' empty tester (lines) or an empty line (testers).
#'
#' @format A data frame with 81 rows and 4 columns:
#' \describe{
#'   \item{Block}{Block (1-3).}
#'   \item{Line}{Line (L1-L6), `NA` for the tester parents.}
#'   \item{Tester}{Tester (T1-T3), `NA` for the line parents.}
#'   \item{Yield_t_ha}{Grain yield, t/ha.}
#' }
#' @examples
#' line_tester(maize_lxt, "Yield_t_ha", "Line", "Tester", block = "Block")
"maize_lxt"

#' North Carolina mating designs
#'
#' Simulated data of the three North Carolina designs, each in two sets
#' with two replicates per set.
#'
#' @format A list of three data frames:
#' \describe{
#'   \item{nc1}{Design I: `Set`, `Rep`, `Male` (three per set), `Female`
#'     (four per male, labelled within the male) and `Yield_t_ha`; 48 rows.}
#'   \item{nc2}{Design II: `Set`, `Rep`, `Male` and `Female` (four of each
#'     per set) and `Yield_t_ha`; 64 rows.}
#'   \item{nc3}{Design III: `Set`, `Rep`, `Male` (six F2 plants per set),
#'     `Tester` (`P1`, `P2`) and `Yield_t_ha`; 48 rows.}
#' }
#' @examples
#' nc_design(nc_designs$nc1, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 1)
"nc_designs"
