# Self-contained HTML report of an analysis: methods paragraph, data, assumptions,
# ANOVA, means with letters, trends, simple effects, figures (inline SVG),
# appendices and citation. The source stays ASCII: symbols are HTML entities.

.esc <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}
# Pre-rounding to 9 decimals makes the same value print the same everywhere (7.1825 computed
# two ways lies on either side of the half and would print 7.182 in one table, 7.183 in another).
.fx <- function(v, d = 3) ifelse(is.na(v) | !is.finite(v), "&mdash;", formatC(round(v, 9), format = "f", digits = d, big.mark = ""))
.fp <- function(p) ifelse(is.na(p), "&mdash;", ifelse(p < 1e-4, "&lt; 0.0001", formatC(p, format = "f", digits = 4)))
.lab_p <- function(p) if (is.na(p)) "" else if (p < 1e-4) "p &lt; 0.0001" else paste0("p = ", formatC(p, format = "f", digits = 4))

# An HTML table; cols: list of list(label, values, num = TRUE/FALSE). Numbered caption.
.html_table <- function(env, cols, caption, row_class = NULL) {
  env$tab <- env$tab + 1
  n <- length(cols[[1]]$values)
  th <- paste0(vapply(cols, function(c) sprintf("<th%s>%s</th>", if (isTRUE(c$num)) " class=\"num\"" else "", c$label), ""), collapse = "")
  rows <- vapply(seq_len(n), function(i) {
    td <- vapply(cols, function(c) sprintf("<td%s>%s</td>", if (isTRUE(c$num)) " class=\"num\"" else "", c$values[i]), "")
    sprintf("<tr%s>%s</tr>", if (!is.null(row_class) && nzchar(row_class[i])) sprintf(" class=\"%s\"", row_class[i]) else "",
            paste0(td, collapse = ""))
  }, "")
  sprintf("<table><caption><b>Table %d.</b> %s</caption><thead><tr>%s</tr></thead><tbody>%s</tbody></table>",
          env$tab, caption, th, paste0(rows, collapse = ""))
}

# ---------- figures embedded in the page ----------

.base64 <- function(raw) {
  chars <- c(LETTERS, letters, 0:9, "+", "/")
  n <- length(raw); pad <- (3 - n %% 3) %% 3
  b <- c(as.integer(raw), rep(0L, pad))
  m <- matrix(b, nrow = 3)
  v <- m[1, ] * 65536L + m[2, ] * 256L + m[3, ]
  idx <- rbind(v %/% 262144L, (v %/% 4096L) %% 64L, (v %/% 64L) %% 64L, v %% 64L) + 1L
  out <- chars[as.vector(idx)]
  if (pad) out[(length(out) - pad + 1):length(out)] <- "="
  paste0(out, collapse = "")
}

# Draws `expr` into an SVG (cairo) or, without cairo, a PNG, and returns the markup to inline.
.inline_figure <- function(env, expr, caption, width = 16, height = 10.5) {
  env$fig <- env$fig + 1
  cairo <- .cairo_ok()
  f <- tempfile(fileext = if (cairo) ".svg" else ".png")
  on.exit(unlink(f))
  ok <- tryCatch({ save_figure(f, expr, width = width, height = height, res = 200); TRUE }, error = function(e) FALSE)
  if (!ok || !file.exists(f)) { env$fig <- env$fig - 1; return("") }
  body <- if (cairo) {
    s <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    s <- sub("^<\\?xml[^>]*>\\s*", "", s)
    # every figure gets its own ids, or the glyphs of one figure would be drawn with those of another
    p <- sprintf("f%d-", env$fig)
    s <- gsub("id=\"", paste0("id=\"", p), s, fixed = TRUE)
    s <- gsub("href=\"#", paste0("href=\"#", p), s, fixed = TRUE)
    s <- gsub("url(#", paste0("url(#", p), s, fixed = TRUE)
    s <- sub("<svg ", "<svg style=\"max-width:100%;height:auto\" ", s, fixed = TRUE)
    sub("width=\"[0-9.]+pt\" height=\"[0-9.]+pt\"", "", s)
  } else {
    sprintf("<img alt=\"Figure %d\" style=\"max-width:100%%\" src=\"data:image/png;base64,%s\">", env$fig,
            .base64(readBin(f, "raw", file.size(f))))
  }
  sprintf("<figure>%s<figcaption><b>Figure %d.</b> %s</figcaption></figure>", body, env$fig, caption)
}

.css <- "
body{font-family:Arial,Helvetica,sans-serif;color:#1b1f2a;max-width:960px;margin:0 auto;padding:28px 32px;line-height:1.5;font-size:14px;background:#fff}
h1{font-size:26px;margin:0 0 4px}h2{font-size:19px;margin:30px 0 8px;border-bottom:2px solid #2f7d4f;padding-bottom:4px;color:#256640}h3{font-size:15px;margin:18px 0 6px}
.meta{color:#5f6b64;margin-bottom:18px}.meta b{color:#1b1f2a}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin:10px 0 16px}caption{caption-side:top;text-align:left;padding:4px 0 6px;font-size:13px}
th,td{padding:5px 8px;border-bottom:1px solid #dde4dc;text-align:left;white-space:nowrap}thead th{border-bottom:2px solid #9aa39d;background:#f4f6f3}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}tr.total td{font-weight:700;border-top:2px solid #9aa39d}tr.dim td{color:#5f6b64}
figure{margin:14px 0 22px;page-break-inside:avoid}figcaption{font-size:12.5px;color:#333;margin-top:6px}
.box{border-left:3px solid #2f7d4f;background:#f4f6f3;padding:10px 14px;margin:10px 0;border-radius:0 6px 6px 0;font-size:13px}.box.warn{border-color:#d98a1c}.box.bad{border-color:#c93a2c}
.grade{display:inline-block;font-weight:800;font-size:22px;padding:2px 12px;border-radius:8px;background:#e9eee7;margin-right:8px}
.small{font-size:12px;color:#5f6b64}code{background:#e9eee7;padding:1px 4px;border-radius:3px;font-size:12px}
pre{background:#f4f6f3;border-radius:8px;padding:10px 12px;font-size:12px;overflow-x:auto}
.footer{margin-top:40px;border-top:1px solid #dde4dc;padding-top:10px;font-size:11.5px;color:#5f6b64}
@media print{body{padding:0;max-width:none;font-size:12px}h2{page-break-after:avoid}table{font-size:11px}a{color:inherit;text-decoration:none}.nobreak{page-break-inside:avoid}}
@page{size:A4;margin:18mm}"

# ---------- content ----------

# "22 October 2026" whatever the locale of the session.
.date_en <- function(d) sprintf("%d %s %s", as.integer(format(d, "%d")), month.name[as.integer(format(d, "%m"))], format(d, "%Y"))

.is_numeric_levels <- function(lv) length(lv) >= 3 && !anyNA(suppressWarnings(as.numeric(lv)))

.pkg_version <- function() tryCatch(as.character(utils::packageVersion("AgriDesign")), error = function(e) "development")

.methods_text <- function(fit, method, alpha, units, with_assumptions) {
  r <- fit$roles; d <- fit$data
  facs <- paste(vapply(r$treatment, function(f) sprintf("%s (%d levels: %s)", .esc(f), nlevels(d[[f]]),
                                                       .esc(paste(levels(d[[f]]), collapse = ", "))), ""), collapse = "; ")
  reps <- if (length(r$block)) sprintf("%d %s", nlevels(d[[r$block[1]]]), if (length(r$block) == 2) "replicates" else "blocks")
          else if (!is.null(r$row)) sprintf("%d rows &times; %d columns", nlevels(d[[r$row]]), nlevels(d[[r$col]]))
          else sprintf("%s replicates per treatment", format(round(fit$n / prod(vapply(r$treatment, function(f) nlevels(d[[f]]), 0)), 1)))
  strata <- any(fit$table$is_error)
  trend <- any(vapply(r$treatment, function(f) .is_numeric_levels(levels(d[[f]])), NA))
  paste0(
    "<p>The experiment was laid out as a <b>", .esc(tolower(substr(fit$design_name, 1, 1))), .esc(substring(sub(" \\(.*", "", fit$design_name), 2)),
    "</b> with ", reps, ". Treatment factor", if (length(r$treatment) > 1) "s were " else " was ", facs,
    ". The response variable was <b>", .esc(fit$response), "</b>", if (nzchar(units)) paste0(" (", .esc(units), ")") else "", ".",
    if (length(r$covariate)) paste0(" ", .esc(paste(r$covariate, collapse = ", ")), if (length(r$covariate) > 1) " were" else " was",
                                    " included as covariate", if (length(r$covariate) > 1) "s", " (analysis of covariance)."),
    if (!is.null(fit$transform)) paste0(" The response was transformed (", .esc(fit$transform$label),
                                        ") before the analysis to meet the assumptions; means are reported back-transformed."),
    if (with_assumptions) paste0(" Residuals were checked for normality (Shapiro&ndash;Wilk), homogeneity of variances (Levene's test",
                                 if (.has_blocking(fit)) ") and additivity (Tukey's one-degree-of-freedom test" else "", ")."),
    " Data were analysed by analysis of variance with the model <code>", .esc(paste(deparse(fit$formula, width.cutoff = 500L), collapse = " ")),
    "</code>, ", if (strata) "testing each term against its own error stratum" else "with a single error term",
    ", and Type ", if (fit$ss_type == 1) "I" else "III", " sums of squares. Treatment means were compared with ",
    .methods[[method]]$name, " at &alpha; = ", alpha,
    if (trend) ", and trends over quantitative factors were examined with orthogonal polynomial contrasts" else "",
    ". All computations were done with the R package AgriDesign version ", .pkg_version(),
    " (Barrera-Guzm&aacute;n, ", format(Sys.Date(), "%Y"), ") in ", .esc(R.version.string), " (R Core Team).</p>")
}

.narrative <- function(fit, method, alpha) {
  t <- fit$table[!fit$table$is_error, ]
  trt <- t[vapply(t$term, function(x) all(strsplit(x, ":", fixed = TRUE)[[1]] %in% fit$roles$treatment), NA), ]
  s <- vapply(seq_len(nrow(trt)), function(i) {
    inter <- grepl(":", trt$term[i], fixed = TRUE)
    sprintf("The %s was %ssignificant (F<sub>%d,%s</sub> = %s, %s).",
            if (inter) paste("interaction", gsub(":", " &times; ", .esc(trt$term[i]), fixed = TRUE)) else paste("effect of", .esc(trt$term[i])),
            if (!is.na(trt$p[i]) && trt$p[i] < alpha) "" else "not ", trt$df[i], format(round(trt$df_error[i], 1)),
            .fx(trt$F[i], 2), .lab_p(trt$p[i]))
  }, "")
  paste0("<p>", paste(c(s, sprintf("The coefficient of variation was %s %%. Means were separated with %s at &alpha; = %s.",
                                   .fx(fit$cv, 1), .methods[[method]]$name, alpha)), collapse = " "), "</p>")
}

#' Write an HTML report of an analysis
#'
#' Builds a self-contained HTML report of a fitted design: a materials and
#' methods paragraph ready to adapt, descriptive statistics, the assumption
#' checks with their grade, the analysis of variance with a draft of the
#' results, the means with their letters, trend analysis of quantitative
#' factors, simple effects of two-factor interactions, figures embedded as
#' vector graphics, appendices with the pairwise comparisons and the data,
#' and how to cite the software. Open it in any web browser and print it to
#' PDF from there (A4 page style included).
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param file Output file (`.html`); `NULL` only returns the HTML.
#' @param title,authors,affiliation,trial,objective Texts for the head of the report.
#' @param units Units of the response, for the methods paragraph.
#' @param method,alpha Mean-separation procedure and significance level (see [compare_means()]).
#' @param control Control level for `method = "dunnett"`.
#' @param assumptions,figures,pairwise,data,cite Include those sections.
#' @param open Open the report in the browser when it is written.
#' @return The HTML text, invisibly.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' f <- tempfile(fileext = ".html")
#' design_report(fit, f, title = "Nitrogen rates in maize", units = "t/ha", figures = FALSE)
#' file.exists(f)
#' @export
design_report <- function(fit, file = NULL, title = NULL, authors = "", affiliation = "", trial = "",
                          objective = "", units = "", method = "tukey", alpha = 0.05, control = NULL,
                          assumptions = TRUE, figures = TRUE, pairwise = TRUE, data = FALSE, cite = TRUE,
                          open = FALSE) {
  .check_fit(fit)
  method <- match.arg(method, names(.methods))
  env <- new.env(); env$tab <- 0; env$fig <- 0
  r <- fit$roles; d <- fit$data
  if (is.null(title)) title <- sprintf("Effect of %s on %s", paste(r$treatment, collapse = " and "), fit$response)
  head <- paste0("<h1>", .esc(title), "</h1><div class=\"meta\">",
                 if (nzchar(authors)) paste0("<b>", .esc(authors), "</b><br>"),
                 if (nzchar(affiliation)) paste0(.esc(affiliation), "<br>"),
                 if (nzchar(trial)) paste0(.esc(trial), "<br>"),
                 "Report generated on ", .date_en(Sys.Date()), " with the R package AgriDesign ", .pkg_version(), "</div>")
  out <- head
  if (nzchar(objective)) out <- c(out, "<h2>Objective</h2>", paste0("<p>", .esc(objective), "</p>"))
  out <- c(out, "<h2>Materials and methods</h2>", .methods_text(fit, method, alpha, units, assumptions))

  # data
  grp <- .treatment_groups(fit); y <- fit$y
  st <- function(v) c(n = length(v), mean = mean(v), sd = if (length(v) > 1) stats::sd(v) else NA,
                      se = if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA, cv = stats::sd(v) / mean(v) * 100,
                      min = min(v), max = max(v))
  ds <- rbind(t(vapply(split(y, grp), st, numeric(7))), All = st(y))
  out <- c(out, "<h2>Data</h2>",
           sprintf("<p>%d observations of %s entered the analysis%s.</p>", fit$n, .esc(fit$response),
                   if (fit$dropped) sprintf(" (%d incomplete rows were excluded)", fit$dropped) else ""),
           .html_table(env, list(
             list(label = paste(.esc(r$treatment), collapse = " &times; "), values = .esc(c(rownames(ds)[-nrow(ds)], "All observations"))),
             list(label = "n", values = ds[, "n"], num = TRUE), list(label = "Mean", values = .fx(ds[, "mean"]), num = TRUE),
             list(label = "SD", values = .fx(ds[, "sd"]), num = TRUE), list(label = "SE", values = .fx(ds[, "se"]), num = TRUE),
             list(label = "CV %", values = .fx(ds[, "cv"], 1), num = TRUE), list(label = "Min", values = .fx(ds[, "min"]), num = TRUE),
             list(label = "Max", values = .fx(ds[, "max"]), num = TRUE)),
             sprintf("Descriptive statistics of %s by treatment.", .esc(fit$response)), c(rep("", nrow(ds) - 1), "total")))

  # assumptions
  if (assumptions) {
    a <- check_assumptions(fit, alpha)
    t <- a$tests[!is.na(a$tests$p) | a$tests$check == "independence", ]
    out <- c(out, "<h2>Assumptions of the analysis of variance</h2>",
             sprintf("<div class=\"box%s\"><span class=\"grade\">%s</span> %s</div>",
                     switch(a$grade, A = , B = "", C = " warn", D = " bad"), a$grade, .esc(.grade_text[[a$grade]])),
             .html_table(env, list(
               list(label = "Test", values = .esc(t$test)),
               list(label = "Statistic", values = paste(.esc(t$statistic_name), "=", .fx(t$statistic))),
               list(label = "df", values = t$df, num = TRUE),
               list(label = "p", values = .fp(t$p), num = TRUE),
               list(label = "Decision", values = ifelse(is.na(t$p), "&mdash;", ifelse(t$p < alpha, "reject H<sub>0</sub>", "do not reject")))),
               sprintf("Assumption checks on the residuals of %s (&alpha; = %s).", .esc(fit$response), alpha)))
    if (nrow(a$suspicious))
      out <- c(out, sprintf("<p class=\"small\">Observations with |studentized residual| &gt; 2.5 (data rows): %s.</p>",
                            paste(sprintf("%s (%s)", a$suspicious$row, .fx(a$suspicious$studentized, 2)), collapse = ", ")))
  }

  # analysis of variance
  tb <- fit$table
  err_lab <- vapply(tb$error, function(e) if (e == "Residuals") "Residual" else .esc(fit$strata[[e]]$label), "")
  out <- c(out, "<h2>Results</h2><h3>Analysis of variance</h3>",
           .html_table(env, list(
             list(label = "Source", values = .esc(c(tb$source, "Residual error", "Total"))),
             list(label = "df", values = c(tb$df, fit$residual[["df"]], fit$total[["df"]]), num = TRUE),
             list(label = "SS", values = .fx(c(tb$ss, fit$residual[["ss"]], fit$total[["ss"]]), 4), num = TRUE),
             list(label = "MS", values = .fx(c(tb$ms, fit$residual[["ms"]], NA), 4), num = TRUE),
             list(label = "F", values = .fx(c(tb$F, NA, NA)), num = TRUE),
             list(label = "p", values = .fp(c(tb$p, NA, NA)), num = TRUE),
             list(label = "", values = c(.stars(tb$p), "", "")),
             list(label = "Tested against", values = c(ifelse(tb$is_error, "", err_lab), "", ""))),
             sprintf("Analysis of variance of %s (%s). CV = %s %%, R<sup>2</sup> = %s, grand mean = %s. Significance: *** p &lt; 0.001, ** p &lt; 0.01, * p &lt; 0.05, . p &lt; 0.1.",
                     .esc(fit$response), .esc(sub(" \\(.*", "", fit$design_name)), .fx(fit$cv, 1), .fx(fit$r_squared), .fx(fit$grand_mean)),
             c(ifelse(tb$is_error, "dim", ""), "dim", "total")),
           .narrative(fit, method, alpha))

  # means
  out <- c(out, "<h3>Treatment means</h3>")
  cmps <- list()
  for (f in r$treatment) {
    cm <- compare_means(fit, f, method = method, alpha = alpha, control = control)
    cmps[[f]] <- cm
    m <- cm$means[order(-cm$means$mean), ]
    cols <- list(list(label = .esc(f), values = .esc(m$label)), list(label = "n", values = m$n, num = TRUE),
                 list(label = "LS mean", values = .fx(m$mean), num = TRUE), list(label = "SE", values = .fx(m$se), num = TRUE))
    if (!is.null(fit$transform)) cols <- c(cols, list(list(label = "Back-transformed", values = .fx(m$back_transformed), num = TRUE)))
    cols <- c(cols, list(list(label = "Group", values = .esc(m$group))))
    msd <- if (length(cm$msd) == 1) sprintf(" Minimum significant difference = %s.", .fx(cm$msd)) else ""
    out <- c(out, .html_table(env, cols, sprintf(
      "Means of %s by %s. %s (%s, &alpha; = %s; error: %s).%s", .esc(fit$response), .esc(f),
      if (method == "dunnett") "* differs from the control" else "Means followed by the same letter are not significantly different",
      cm$method_name, alpha, .esc(cm$error$label), msd)))
    if (.is_numeric_levels(levels(d[[f]]))) {
      tr <- trend(fit, f, degree = 3, alpha = alpha)
      cp <- tr$components
      out <- c(out, .html_table(env, list(
        list(label = "Component", values = cp$component), list(label = "Coefficients", values = cp$coefficients),
        list(label = "SS", values = .fx(cp$ss, 4), num = TRUE), list(label = "F", values = .fx(cp$F), num = TRUE),
        list(label = "p", values = .fp(cp$p), num = TRUE), list(label = "", values = .stars(cp$p))),
        sprintf("Orthogonal polynomial contrasts for %s.%s", .esc(f),
                if (!is.null(tr$optimum)) sprintf(" The fitted quadratic has its %s at %s = %s.", attr(tr$optimum, "kind"), .esc(f), .fx(tr$optimum[["x"]], 1)) else "")))
    }
  }
  pairs2 <- if (length(r$treatment) >= 2) utils::combn(r$treatment, 2, simplify = FALSE) else list()
  # simple effects in both directions: B within each level of A, and A within each level of B
  slices <- unlist(lapply(pairs2, function(pr) list(pr, rev(pr))), recursive = FALSE)
  for (pr in slices) {
    if (!any(vapply(fit$terms, function(t) setequal(t$factors, pr), NA))) next
    se <- compare_means(fit, pr[2], by = pr[1], method = if (method == "dunnett") "tukey" else method, alpha = alpha)
    m <- se$means
    la <- unique(m[[pr[1]]]); lb <- unique(m$label)
    cols <- c(list(list(label = paste0(.esc(pr[1]), " \\ ", .esc(pr[2])), values = .esc(la))),
              lapply(lb, function(b) list(label = .esc(b), num = TRUE, values = vapply(la, function(a) {
                k <- m[m[[pr[1]]] == a & m$label == b, ]; paste(.fx(k$mean), .esc(k$group)) }, ""))))
    out <- c(out, .html_table(env, cols, sprintf(
      "Cell means of %s for %s &times; %s. Letters compare the levels of %s within each level of %s (%s, &alpha; = %s; error: %s).",
      .esc(fit$response), .esc(pr[1]), .esc(pr[2]), .esc(pr[2]), .esc(pr[1]), se$method_name, alpha, .esc(se$error$label))))
  }

  # figures
  if (figures) {
    figs <- character(0)
    for (f in r$treatment) {
      cm <- cmps[[f]]
      figs <- c(figs, .inline_figure(env, plot(cm), sprintf("Means of %s by %s (&plusmn; SE) with the groups of %s.", .esc(fit$response), .esc(f), cm$method_name)))
      if (.is_numeric_levels(levels(d[[f]])))
        figs <- c(figs, .inline_figure(env, plot(trend(fit, f, alpha = alpha)), sprintf("Response of %s to %s: means &plusmn; SE and fitted polynomial.", .esc(fit$response), .esc(f))))
    }
    for (pr in pairs2) if (any(vapply(fit$terms, function(t) setequal(t$factors, pr), NA)))
      figs <- c(figs, .inline_figure(env, plot_interaction(fit, pr[2], pr[1]), sprintf("Interaction %s &times; %s: cell means &plusmn; SE.", .esc(pr[1]), .esc(pr[2]))))
    if (any(fit$table$is_error)) figs <- c(figs, .inline_figure(env, plot_ss_partition(fit), "Partition of the total sum of squares.", height = 8))
    if (assumptions) figs <- c(figs, .inline_figure(env, plot(fit), "Residual diagnostics of the model.", width = 24, height = 15))
    if (!is.null(r$row) && !is.null(r$col))
      figs <- c(figs, .inline_figure(env, plot_field(fit, value = "residuals"), "Residuals on the field map (rows &times; columns).", width = 14))
    figs <- figs[nzchar(figs)]
    if (length(figs)) out <- c(out, sprintf("<div class=\"nobreak\"><h2>Figures</h2>%s</div>", figs[1]), figs[-1])
  }

  if (pairwise) {
    out <- c(out, "<h2>Appendix A &middot; Pairwise comparisons</h2>")
    for (f in r$treatment) {
      p <- cmps[[f]]$pairs
      out <- c(out, .html_table(env, list(
        list(label = "Level i", values = .esc(p$level_i)), list(label = "Level j", values = .esc(p$level_j)),
        list(label = "Difference", values = .fx(p$diff), num = TRUE), list(label = "SE", values = .fx(p$se), num = TRUE),
        list(label = "Statistic", values = .fx(p$statistic), num = TRUE), list(label = "adj. p", values = .fp(p$p_adj), num = TRUE),
        list(label = "", values = ifelse(p$different, "*", "ns"))),
        sprintf("Pairwise comparisons for %s (%s).", .esc(f), cmps[[f]]$method_name)))
    }
  }
  if (data) {
    dd <- d; dd[[fit$response]] <- fit$y
    out <- c(out, "<h2>Appendix B &middot; Data</h2>",
             .html_table(env, lapply(names(dd), function(n) list(label = .esc(n), values = .esc(format(dd[[n]])), num = is.numeric(dd[[n]]))),
                         "Data used in the analysis."))
  }
  if (cite) {
    yr <- format(Sys.Date(), "%Y")
    out <- c(out, "<h2>How to cite</h2><p>If this analysis is published, please cite the software:</p>",
             sprintf("<p style=\"padding-left:2em;text-indent:-2em\">Barrera-Guzm&aacute;n, L. &Aacute;. (%s). <i>AgriDesign: Design and Analysis of Agricultural Experiments</i> (R package version %s). <a href=\"https://github.com/luisangelbg/AgriDesign-R\">https://github.com/luisangelbg/AgriDesign-R</a></p>", yr, .pkg_version()),
             sprintf("<p style=\"padding-left:2em;text-indent:-2em\">R Core Team (%s). <i>R: A Language and Environment for Statistical Computing</i>. R Foundation for Statistical Computing, Vienna, Austria. <a href=\"https://www.R-project.org/\">https://www.R-project.org/</a></p>", R.version$year))
  }
  out <- c(out, "<div class=\"footer\">AgriDesign &middot; design and analysis of agricultural experiments. Figures are embedded as vector graphics and can be extracted from this file.</div>")
  html <- paste0("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>",
                 .esc(title), "</title><style>", .css, "</style></head><body>", paste(out, collapse = "\n"), "</body></html>")
  if (!is.null(file)) {
    con <- file(file, open = "w", encoding = "UTF-8"); on.exit(close(con))
    writeLines(html, con, useBytes = FALSE)
    if (open) utils::browseURL(normalizePath(file))
  }
  attr(html, "tables") <- env$tab; attr(html, "figures") <- env$fig
  invisible(html)
}
