# Publication figures in base graphics: residual diagnostics, Box-Cox profile,
# means with letters, simple effects, interaction, response curve, rank-based
# groups, differences with intervals, partition of the sums of squares, cell
# means and field maps; save_figure() writes them at a given size and resolution.

.agri_palette <- c("#2f7d4f", "#c8842a", "#2b7bb9", "#b5432f", "#7a5195",
                   "#5aa469", "#e0b04c", "#3b9db3", "#8c5a3c", "#6d6e71")

.ramps <- list(
  ylgn = c("#fcfacb", "#bde07e", "#47944c", "#034741"),
  rdbu = c("#023869", "#97c5de", "#f6efe6", "#e5a692", "#902822"))

#' Colours of the AgriDesign figures
#'
#' The categorical palette used by the figures of the package (forest green,
#' ochre, blue, brick red, plum...), recycled to `n` colours.
#'
#' @param n Number of colours.
#' @return A character vector of colours.
#' @examples
#' agri_colors(4)
#' @export
agri_colors <- function(n) rep_len(.agri_palette, n)

.ramp <- function(name, n = 100) grDevices::colorRampPalette(.ramps[[name]])(n)

.darken <- function(col, f = 0.25) {
  m <- grDevices::col2rgb(col) / 255 * (1 - f)
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

.lighten <- function(col, f = 0.6) {
  m <- grDevices::col2rgb(col) / 255
  m <- m + (1 - m) * f
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

# Text colour that reads on top of a fill.
.on_color <- function(col) {
  m <- grDevices::col2rgb(col) / 255
  ifelse(0.299 * m[1, ] + 0.587 * m[2, ] + 0.114 * m[3, ] > 0.55, "#1b1b1b", "#ffffff")
}

.style_par <- function(mar = c(4.3, 4.6, 2.8, 1.2), ...) {
  graphics::par(mar = mar, mgp = c(2.7, 0.7, 0), tcl = -0.3, las = 1, cex.axis = 0.9, bty = "l", ...)
}

.note <- function(text) graphics::mtext(text, side = 3, adj = 1, line = 0.2, cex = 0.72, col = "#5f6368")

.long_labels <- function(labels, k) sum(nchar(labels)) > 70 || max(nchar(labels)) > 14 && k > 4

# Error half-width of each mean.
.err_width <- function(m, error, df) {
  switch(error, se = m$se, ci = stats::qt(0.975, df) * m$se, none = rep(0, nrow(m)))
}

# Plus-minus, middle dot and superscripts, built from their code points so that
# the source stays ASCII.
.pm <- intToUtf8(177); .dot <- intToUtf8(183); .sup <- c(intToUtf8(178), intToUtf8(179), "^4", "^5")
.sep <- paste0("  ", .dot, "  ")
.pch_series <- c(21, 22, 24, 23, 25)

# ", p = 0.012" or ", p < 2e-16"
.p_text <- function(p) {
  f <- format.pval(p, digits = 3)
  ifelse(substr(f, 1, 1) == "<", paste0(", p ", f), paste0(", p = ", f))
}

.error_label <- function(error) switch(error, se = paste("Mean", .pm, "SE"), ci = paste("Mean", .pm, "95 % CI"), none = "Mean")

.ylab_of <- function(x) {
  if (is.null(x$transform)) x$response else paste0(x$response, " (", x$transform$label, ")")
}

# ---------- means with letters ----------

#' Plot treatment means with their letters
#'
#' Draws the means of a comparison from [compare_means()] as bars or points
#' with error bars, or the observations as box plots or points, and writes
#' the letters of the groups above them. With simple effects (`by`), the
#' means are grouped by the levels of the other factor and the letters are
#' valid within each group.
#'
#' @param x An `"ad_comparison"` object.
#' @param style `"bars"`, `"points"`, `"box"` (box plots of the
#'   observations) or `"strip"` (every observation, with the mean and its
#'   error bar).
#' @param error Error bars: `"se"` (standard error of the mean), `"ci"`
#'   (95 % confidence interval with the error degrees of freedom) or `"none"`.
#' @param letters Write the letters of the groups.
#' @param values Write the value of each mean.
#' @param sort Sort the levels by their mean, highest first.
#' @param col Fill colours, one per level (default: [agri_colors()]).
#' @param main,xlab,ylab Titles; `NULL` builds them from the analysis.
#' @param ylim Range of the vertical axis.
#' @param letter_cex Size of the letters.
#' @param digits Decimals of the values written with `values = TRUE`.
#' @param ... Further graphical parameters passed to [graphics::plot.window()].
#' @return Invisibly, a data frame with the position, mean, error-bar limits
#'   and letters of every level.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' cmp <- compare_means(fit, "Nitrogen_kg_ha")
#' plot(cmp)
#' plot(cmp, style = "box")
#'
#' fit <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
#'                     block = "Block", design = "split_rcbd")
#' plot(compare_means(fit, "Variety", by = "Irrigation"))
#' @export
plot.ad_comparison <- function(x, style = c("bars", "points", "box", "strip"), error = c("se", "ci", "none"),
                               letters = TRUE, values = FALSE, sort = FALSE, col = NULL,
                               main = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
                               letter_cex = 1.1, digits = 2, ...) {
  style <- match.arg(style); error <- match.arg(error)
  if (!is.null(x$by)) return(invisible(.plot_grouped(x, error, letters, col, main, xlab, ylab, ylim, letter_cex, ...)))
  m <- x$means
  if (sort) m <- m[order(-m$mean), ]
  k <- nrow(m)
  col <- if (is.null(col)) agri_colors(k) else rep_len(col, k)
  e <- .err_width(m, error, x$error$df)
  lo <- m$mean - e; hi <- m$mean + e
  obs <- split(x$obs$y, factor(x$obs$label, levels = m$label))
  top <- if (style %in% c("box", "strip")) pmax(vapply(obs, max, 0), hi) else hi
  bottom <- if (style %in% c("box", "strip")) pmin(vapply(obs, min, 0), lo) else lo
  if (is.null(ylim)) {
    r <- range(c(bottom, top)); pad <- diff(r) * 0.14
    ylim <- if (style == "bars") c(min(0, r[1]), r[2] + pad + if (values) pad * 0.6 else 0)
            else c(r[1] - diff(r) * 0.05, r[2] + pad + if (values) pad * 0.6 else 0)
  }
  long <- .long_labels(m$label, k)
  op <- .style_par(mar = c(if (long) 7.5 else 4.3, 4.6, 2.8, 1.2)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.4, k + 0.6), ylim = ylim, xaxs = "i", ...)
  graphics::abline(h = pretty(ylim), col = "#ececec", lwd = 0.8)
  xs <- seq_len(k); w <- 0.32
  if (style == "bars") {
    base <- max(ylim[1], 0)
    graphics::rect(xs - w, base, xs + w, m$mean, col = col, border = .darken(col))
  } else if (style == "box") {
    graphics::boxplot(obs, at = xs, add = TRUE, axes = FALSE, col = .lighten(col, 0.55), border = .darken(col, 0.35),
                      boxwex = 0.55, outpch = 21, outbg = col, staplewex = 0.4)
  } else if (style == "strip") {
    for (i in xs) {
      v <- obs[[i]]; n <- length(v)
      off <- if (n > 1) (seq_len(n) - (n + 1) / 2) / n * 0.32 else 0
      graphics::points(i - 0.08 + off, v, pch = 21, bg = grDevices::adjustcolor(col[i], 0.55), col = .darken(col[i]), cex = 1.05)
    }
  }
  # in the strip style the mean and its error bar sit beside the observations
  mx <- if (style == "strip") xs + 0.26 else xs
  if (error != "none" && style != "box") {
    graphics::arrows(mx, lo, mx, hi, angle = 90, code = 3, length = 0.045, lwd = 1.4, col = "#1b1b1b")
  }
  if (style %in% c("points", "strip")) {
    graphics::points(mx, m$mean, pch = if (style == "strip") 23 else 21, bg = if (style == "strip") "#1b1b1b" else col,
                     col = "#ffffff", cex = if (style == "strip") 1.2 else 1.7, lwd = 1.2)
  }
  ly <- top + diff(ylim) * 0.035
  if (values) {
    graphics::text(xs, ly, formatC(m$mean, format = "f", digits = digits), adj = c(0.5, 0), cex = 0.8, col = "#3c4043")
    ly <- ly + diff(ylim) * 0.055
  }
  if (letters && !is.null(m$group)) {
    lab <- if (x$method == "dunnett") ifelse(m$group == "ns", "", ifelse(m$group == "control", "C", m$group)) else m$group
    graphics::text(xs, ly, lab, adj = c(0.5, 0), font = 2, cex = letter_cex)
  }
  graphics::axis(2); graphics::box()
  graphics::axis(1, at = xs, labels = m$label, tick = FALSE, las = if (long) 2 else 1, line = -0.4)
  graphics::title(main = if (is.null(main)) paste(x$response, "by", paste(x$factor, collapse = " x ")) else main,
                  xlab = if (is.null(xlab)) (if (long) "" else paste(x$factor, collapse = " x ")) else xlab,
                  ylab = if (is.null(ylab)) .ylab_of(x) else ylab, font.main = 1, cex.main = 1.05)
  .note(paste0(if (style %in% c("box", "strip")) "Observations" else .error_label(error),
               if (letters) paste0(.sep, "letters: ", x$method_name) else ""))
  invisible(data.frame(label = m$label, x = xs, mean = m$mean, lower = lo, upper = hi,
                       letters = if (is.null(m$group)) NA else m$group, letter_y = ly, stringsAsFactors = FALSE))
}

.plot_grouped <- function(x, error, letters, col, main, xlab, ylab, ylim, letter_cex, ...) {
  m <- x$means
  g_lv <- unique(m[[x$by]]); f_lv <- unique(m$label)
  nb <- length(g_lv); nf <- length(f_lv)
  col <- if (is.null(col)) agri_colors(nf) else rep_len(col, nf)
  e <- .err_width(m, error, x$error$df)
  gap <- 1; w <- 0.8 / nf
  m$x <- (match(m[[x$by]], g_lv) - 1) * (1 + gap * 0.25) + (match(m$label, f_lv) - (nf + 1) / 2) * w + 1
  hi <- m$mean + e
  if (is.null(ylim)) { r <- range(c(0, hi)); ylim <- c(r[1], r[2] + diff(r) * 0.14 + diff(r) * 0.08) }
  op <- .style_par(); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = range(m$x) + c(-w * 1.2, w * 1.2), ylim = ylim, xaxs = "i", ...)
  graphics::abline(h = pretty(ylim), col = "#ececec", lwd = 0.8)
  fc <- col[match(m$label, f_lv)]
  graphics::rect(m$x - w * 0.45, max(ylim[1], 0), m$x + w * 0.45, m$mean, col = fc, border = .darken(fc))
  if (error != "none") graphics::arrows(m$x, m$mean - e, m$x, hi, angle = 90, code = 3, length = 0.035, lwd = 1.2)
  ly <- hi + diff(ylim) * 0.03
  if (letters) graphics::text(m$x, ly, m$group, adj = c(0.5, 0), font = 2, cex = letter_cex * 0.9)
  centers <- tapply(m$x, factor(m[[x$by]], levels = g_lv), mean)
  graphics::axis(2); graphics::box()
  graphics::axis(1, at = centers, labels = g_lv, tick = FALSE, line = -0.4)
  graphics::legend("topleft", legend = f_lv, fill = col, border = .darken(col), bty = "n", cex = 0.85,
                   ncol = min(nf, 4), title = x$factor, inset = c(0.01, 0))
  graphics::title(main = if (is.null(main)) paste0(x$response, ": ", x$factor, " within each ", x$by) else main,
                  xlab = if (is.null(xlab)) x$by else xlab, ylab = if (is.null(ylab)) .ylab_of(x) else ylab,
                  font.main = 1, cex.main = 1.05)
  .note(paste0(.error_label(error), .sep, "letters within each ", x$by, ": ", x$method_name))
  m$letter_y <- ly
  m
}

# ---------- interaction ----------

#' Interaction plot of two treatment factors
#'
#' Least-squares means of every combination of `x` (horizontal axis) and
#' `trace` (one line per level), with error bars.
#'
#' @inheritParams plot.ad_comparison
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param x,trace Treatment factors on the horizontal axis and as lines.
#' @return Invisibly, the cell means plotted.
#' @examples
#' fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
#'                     block = "Block")
#' plot_interaction(fit, "Fertilizer", "Variety")
#' @export
plot_interaction <- function(fit, x, trace, error = c("se", "ci", "none"), col = NULL,
                             main = NULL, xlab = NULL, ylab = NULL, ylim = NULL) {
  .check_fit(fit); .check_factors(fit, c(x, trace))
  error <- match.arg(error)
  cells <- ls_means(fit, c(x, trace))
  xl <- levels(fit$data[[x]]); tl <- levels(fit$data[[trace]])
  col <- if (is.null(col)) agri_colors(length(tl)) else rep_len(col, length(tl))
  df <- fit$strata[[.error_of(fit, c(x, trace))]]$df
  if (is.null(df)) df <- fit$residual[["df"]]
  e <- .err_width(cells, error, df)
  cells$px <- match(cells[[x]], xl) + (match(cells[[trace]], tl) - (length(tl) + 1) / 2) * 0.04
  cells$err <- e
  if (is.null(ylim)) { r <- range(c(cells$mean - e, cells$mean + e)); ylim <- r + c(-0.05, 0.2) * diff(r) }
  op <- .style_par(); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.7, length(xl) + 0.3), ylim = ylim)
  graphics::abline(h = pretty(ylim), col = "#ececec", lwd = 0.8)
  for (j in seq_along(tl)) {
    s <- cells[cells[[trace]] == tl[j], ]; s <- s[order(match(s[[x]], xl)), ]
    if (error != "none") graphics::arrows(s$px, s$mean - s$err, s$px, s$mean + s$err, angle = 90, code = 3,
                                          length = 0.035, col = col[j], lwd = 1.2)
    graphics::lines(s$px, s$mean, col = col[j], lwd = 2)
    graphics::points(s$px, s$mean, pch = .pch_series[(j - 1) %% 5 + 1], bg = col[j], col = "#ffffff", cex = 1.5, lwd = 1)
  }
  graphics::axis(1, at = seq_along(xl), labels = xl); graphics::axis(2); graphics::box()
  graphics::legend("topleft", legend = tl, col = col, pt.bg = col, pch = .pch_series[(seq_along(tl) - 1) %% 5 + 1],
                   lwd = 2, bty = "n", cex = 0.85, title = trace, ncol = min(length(tl), 4))
  graphics::title(main = if (is.null(main)) paste0(fit$response, ": ", x, " x ", trace) else main,
                  xlab = if (is.null(xlab)) x else xlab, ylab = if (is.null(ylab)) .ylab_of(fit) else ylab,
                  font.main = 1, cex.main = 1.05)
  p <- fit$table$p[fit$table$term %in% c(paste(x, trace, sep = ":"), paste(trace, x, sep = ":"))]
  .note(paste0(.error_label(error), if (length(p)) paste0(.sep, "interaction p = ", format.pval(p, digits = 3)) else ""))
  invisible(cells)
}

#' Heat map of cell means
#'
#' Least-squares means of every combination of two treatment factors as a
#' coloured table with the values written in the cells.
#'
#' @inheritParams plot_interaction
#' @param rows,cols Factors in the rows and columns of the map.
#' @param digits Decimals of the values.
#' @param main Title; `NULL` builds it from the analysis.
#' @return Invisibly, the cell means plotted.
#' @examples
#' fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
#'                     block = "Block")
#' plot_cell_means(fit, "Variety", "Fertilizer")
#' @export
plot_cell_means <- function(fit, rows, cols, digits = 2, main = NULL) {
  .check_fit(fit); .check_factors(fit, c(rows, cols))
  cells <- ls_means(fit, c(rows, cols))
  rl <- levels(fit$data[[rows]]); cl <- levels(fit$data[[cols]])
  v <- cells$mean
  pal <- .ramp("ylgn")
  fill <- pal[1 + round((v - min(v)) / max(diff(range(v)), 1e-12) * 99)]
  op <- .style_par(mar = c(4.3, 6.5, 2.8, 5)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, length(cl) + 0.5), ylim = c(length(rl) + 0.5, 0.5), xaxs = "i", yaxs = "i")
  ci <- match(cells[[cols]], cl); ri <- match(cells[[rows]], rl)
  graphics::rect(ci - 0.48, ri - 0.48, ci + 0.48, ri + 0.48, col = fill, border = NA)
  graphics::text(ci, ri, formatC(v, format = "f", digits = digits), col = .on_color(fill), font = 2, cex = 0.95)
  graphics::axis(1, at = seq_along(cl), labels = cl, tick = FALSE)
  graphics::axis(2, at = seq_along(rl), labels = rl, tick = FALSE)
  .colorbar(pal, range(v), digits)
  graphics::title(main = if (is.null(main)) paste0("Mean ", fit$response, " by ", rows, " and ", cols) else main,
                  xlab = cols, ylab = "", font.main = 1, cex.main = 1.05)
  graphics::mtext(rows, side = 2, line = 5, las = 0)
  invisible(cells)
}

# Vertical colour bar in the right margin of the current plot.
.colorbar <- function(pal, rng, digits) {
  u <- graphics::par("usr"); fx <- diff(u[1:2])
  x0 <- u[2] + fx * 0.04; x1 <- x0 + fx * 0.035
  yb <- seq(min(u[3:4]), max(u[3:4]), length.out = length(pal) + 1)
  if (u[3] > u[4]) pal <- rev(pal)
  graphics::rect(x0, yb[-length(yb)], x1, yb[-1], col = pal, border = NA, xpd = NA)
  at <- c(0, 0.5, 1)
  ys <- min(u[3:4]) + at * abs(diff(u[3:4]))
  vals <- if (u[3] > u[4]) rev(rng[1] + at * diff(rng)) else rng[1] + at * diff(rng)
  graphics::text(x1 + fx * 0.012, ys, formatC(vals, format = "f", digits = digits), adj = 0, cex = 0.75, xpd = NA)
}

# ---------- response curve ----------

#' @export
plot.ad_trend <- function(x, error = c("se", "ci", "none"), main = NULL, xlab = NULL, ylab = NULL, col = "#2f7d4f", ...) {
  error <- match.arg(error)
  m <- x$means
  e <- .err_width(m, error, x$error$df)
  xs <- seq(min(m$x), max(m$x), length.out = 200)
  ys <- drop(outer(xs, seq_along(x$coef) - 1, `^`) %*% x$coef)
  r <- range(c(m$mean - e, m$mean + e, ys))
  op <- .style_par(); on.exit(graphics::par(op))
  graphics::plot(m$x, m$mean, type = "n", ylim = r + c(-0.05, 0.12) * diff(r), xlab = "", ylab = "", ...)
  graphics::abline(h = pretty(r), col = "#ececec", lwd = 0.8)
  graphics::lines(xs, ys, col = col, lwd = 2.2)
  if (error != "none") graphics::arrows(m$x, m$mean - e, m$x, m$mean + e, angle = 90, code = 3, length = 0.04, lwd = 1.3)
  graphics::points(m$x, m$mean, pch = 21, bg = col, col = "#ffffff", cex = 1.7, lwd = 1.2)
  if (!is.null(x$optimum)) {
    o <- x$optimum
    if (o[["x"]] >= min(m$x) && o[["x"]] <= max(m$x)) {
      graphics::segments(o[["x"]], graphics::par("usr")[3], o[["x"]], o[["y"]], lty = 2, col = "#5f6368")
      graphics::points(o[["x"]], o[["y"]], pch = 4, cex = 1.4, lwd = 2, col = "#b5432f")
    }
  }
  graphics::legend("topleft", legend = .equation(x$coef, x$r_squared), bty = "n", cex = 0.85)
  graphics::title(main = if (is.null(main)) paste("Response of", x$response, "to", x$factor) else main,
                  xlab = if (is.null(xlab)) x$factor else xlab, ylab = if (is.null(ylab)) x$response else ylab,
                  font.main = 1, cex.main = 1.05)
  .note(paste0(.error_label(error), if (!is.null(x$optimum)) paste0(.sep, sprintf("%s at x = %.4g", attr(x$optimum, "kind"), x$optimum[["x"]])) else ""))
  invisible(m)
}

.equation <- function(b, r2) {
  s <- format(signif(b[1], 4))
  for (i in seq_along(b)[-1]) s <- paste0(s, if (b[i] < 0) " - " else " + ", format(signif(abs(b[i]), 4)),
                                         " x", if (i > 2) .sup[i - 2] else "")
  sprintf(paste0("y = %s   R", .sup[1], " = %.3f"), s, r2)
}

# ---------- rank-based groups ----------

#' @export
plot.ad_nonparametric <- function(x, col = NULL, main = NULL, xlab = NULL, ylab = NULL, letter_cex = 1.1, ...) {
  g <- x$groups
  k <- nrow(g); col <- if (is.null(col)) agri_colors(k) else rep_len(col, k)
  r <- range(c(0, g$median)); ylim <- c(r[1], r[2] + diff(r) * 0.15)
  long <- .long_labels(g$level, k)
  op <- .style_par(mar = c(if (long) 7.5 else 4.3, 4.6, 2.8, 1.2)); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(xlim = c(0.4, k + 0.6), ylim = ylim, xaxs = "i", ...)
  graphics::abline(h = pretty(ylim), col = "#ececec", lwd = 0.8)
  xs <- seq_len(k)
  graphics::rect(xs - 0.32, 0, xs + 0.32, g$median, col = col, border = .darken(col))
  if (!is.null(g$group)) graphics::text(xs, g$median + diff(ylim) * 0.03, g$group, adj = c(0.5, 0), font = 2, cex = letter_cex)
  graphics::axis(2); graphics::box()
  graphics::axis(1, at = xs, labels = g$level, tick = FALSE, las = if (long) 2 else 1, line = -0.4)
  graphics::title(main = if (is.null(main)) paste("Median", x$response, "by", x$factor) else main,
                  xlab = if (is.null(xlab)) (if (long) "" else x$factor) else xlab,
                  ylab = if (is.null(ylab)) paste("Median", x$response) else ylab, font.main = 1, cex.main = 1.05)
  .note(paste0(sub(" \\(.*", "", x$test$name), if (nzchar(x$posthoc)) paste0(.sep, "groups: ", x$posthoc) else ""))
  invisible(g)
}

# ---------- differences ----------

#' Plot the pairwise differences with their intervals
#'
#' Every pairwise difference of a comparison with its interval: the critical
#' difference of the procedure (simultaneous for Tukey, Scheffe, Dunnett...)
#' or, for Holm and Games-Howell, the 95 % t interval. Differences declared
#' significant are drawn in red.
#'
#' @param x An `"ad_comparison"` object from [compare_means()].
#' @param sort Sort by the size of the difference.
#' @param main,xlab Titles.
#' @param col Colours of non-significant and significant differences.
#' @return Invisibly, the differences and interval limits plotted.
#' @examples
#' fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
#' plot_differences(compare_means(fit, "Variety"))
#' @export
plot_differences <- function(x, sort = TRUE, main = NULL, xlab = "Difference between means",
                             col = c("#2f7d4f", "#b5432f")) {
  if (!inherits(x, "ad_comparison")) stop("'x' must come from compare_means()", call. = FALSE)
  p <- x$pairs
  hw <- ifelse(is.na(p$critical), stats::qt(0.975, x$error$df) * p$se, p$critical)
  p$lower <- p$diff - hw; p$upper <- p$diff + hw
  p$label <- paste(p$level_i, "-", p$level_j, if (!is.null(x$by)) paste0("  (", p[[x$by]], ")") else "")
  if (sort) p <- p[order(p$diff), ]
  n <- nrow(p)
  op <- .style_par(mar = c(4.3, min(22, 1 + max(nchar(p$label)) * 0.55), 2.8, 4)); on.exit(graphics::par(op))
  xl <- range(c(0, p$lower, p$upper))
  graphics::plot.new(); graphics::plot.window(xlim = xl, ylim = c(0.5, n + 0.5))
  graphics::abline(v = pretty(xl), col = "#ececec", lwd = 0.8)
  graphics::abline(v = 0, lty = 2, col = "#3c4043")
  cc <- ifelse(p$different, col[2], col[1])
  graphics::arrows(p$lower, seq_len(n), p$upper, seq_len(n), angle = 90, code = 3, length = 0.04, col = cc, lwd = 2)
  graphics::points(p$diff, seq_len(n), pch = 19, col = cc, cex = 1.1)
  graphics::axis(1); graphics::axis(2, at = seq_len(n), labels = p$label, tick = FALSE, cex.axis = 0.82)
  graphics::box()
  graphics::title(main = if (is.null(main)) paste("Differences of", x$response) else main, xlab = xlab, font.main = 1, cex.main = 1.05)
  .note(paste0(x$method_name, ", alpha = ", x$alpha, .sep, "red: different"))
  invisible(p)
}

# ---------- partition of the sums of squares ----------

#' Partition of the total sum of squares
#'
#' A bar split into the share of the total sum of squares taken by every
#' source of the analysis of variance, including the error strata and the
#' residual.
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param col Colours of the sources.
#' @param main Title.
#' @return Invisibly, a data frame with the sum of squares and percentage of
#'   every source.
#' @examples
#' fit <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
#'                     block = "Block", design = "split_rcbd")
#' plot_ss_partition(fit)
#' @export
plot_ss_partition <- function(fit, col = NULL, main = NULL) {
  .check_fit(fit)
  t <- fit$table
  rows <- data.frame(source = c(t$source, "Residual"), ss = c(t$ss, fit$residual[["ss"]]),
                     df = c(t$df, fit$residual[["df"]]), p = c(t$p, NA), stringsAsFactors = FALSE)
  rows$pct <- rows$ss / sum(rows$ss) * 100
  n <- nrow(rows); col <- if (is.null(col)) agri_colors(n) else rep_len(col, n)
  op <- .style_par(mar = c(1, 1, 2.8, 1)); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(xlim = c(0, 100), ylim = c(-n * 1.1 - 1, 3.2), xaxs = "i")
  x0 <- c(0, cumsum(rows$pct))[seq_len(n)]
  graphics::rect(x0, 1, x0 + rows$pct, 3, col = col, border = "#ffffff")
  big <- rows$pct > 5
  graphics::text((x0 + rows$pct / 2)[big], 2, sprintf("%.1f %%", rows$pct[big]), col = .on_color(col[big]), font = 2, cex = 0.85)
  ly <- -seq_len(n) * 1.1 + 0.2
  graphics::rect(0, ly - 0.35, 2.2, ly + 0.35, col = col, border = NA)
  graphics::text(3.2, ly, sprintf("%s: SS = %s (%.1f %%), df = %s%s", rows$source, format(rows$ss, digits = 4), rows$pct,
                                  format(rows$df, digits = 3),
                                  ifelse(is.na(rows$p), "", .p_text(rows$p))),
                 adj = 0, cex = 0.85)
  graphics::title(main = if (is.null(main)) paste("Partition of the total SS of", fit$response) else main, font.main = 1, cex.main = 1.05)
  invisible(rows)
}

# ---------- field map ----------

#' Field map of the response or of the residuals
#'
#' Draws the plots of the trial on their field rows and columns, coloured by
#' the response or by the residuals of the model (a diverging scale centred
#' at zero), and labelled with the treatment. Spatial patterns of the
#' residuals point to fertility gradients that blocking did not remove.
#'
#' @param fit An `"ad_anova"` object from [design_anova()].
#' @param row,col Names of the field row and column columns. By default the
#'   row and column of a Latin square.
#' @param data The data frame given to [design_anova()], needed when the
#'   field row and column were not roles of the design.
#' @param value `"response"` or `"residuals"`.
#' @param labels Write the treatment of each plot.
#' @param values Write the value of each plot.
#' @param digits Decimals of the values.
#' @param main Title.
#' @return Invisibly, a data frame with the row, column, treatment and value
#'   of every plot.
#' @examples
#' fit <- design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column")
#' plot_field(fit, value = "residuals")
#' @export
plot_field <- function(fit, row = NULL, col = NULL, data = NULL, value = c("response", "residuals"),
                       labels = TRUE, values = TRUE, digits = 2, main = NULL) {
  .check_fit(fit); value <- match.arg(value)
  row <- if (is.null(row)) fit$roles$row else row
  col <- if (is.null(col)) fit$roles$col else col
  if (is.null(row) || is.null(col)) stop("give the field 'row' and 'col' columns", call. = FALSE)
  src <- if (all(c(row, col) %in% names(fit$data))) fit$data else {
    if (is.null(data)) stop("pass 'data': the field row and column are not roles of the design", call. = FALSE)
    .check_columns(data, c(row, col))
    as.data.frame(data)[as.integer(rownames(fit$data)), c(row, col)]
  }
  rr <- as.character(src[[row]]); cc <- as.character(src[[col]])
  rl <- .levels_of(src[[row]]); cl <- .levels_of(src[[col]])
  v <- if (value == "response") fit$y else fit$qr_fit$residuals
  trt <- if (length(fit$roles$treatment) == 1) as.character(fit$data[[fit$roles$treatment]])
         else do.call(paste, c(lapply(fit$data[fit$roles$treatment], as.character), sep = "\n"))
  if (value == "residuals") {
    pal <- .ramp("rdbu"); m <- max(abs(v)); rng <- c(-m, m)
  } else { pal <- .ramp("ylgn"); rng <- range(v) }
  fill <- pal[1 + round((v - rng[1]) / max(diff(rng), 1e-12) * 99)]
  op <- .style_par(mar = c(4.3, 4.6, 2.8, 5)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, length(cl) + 0.5), ylim = c(length(rl) + 0.5, 0.5), xaxs = "i", yaxs = "i")
  xi <- match(cc, cl); yi <- match(rr, rl)
  graphics::rect(xi - 0.47, yi - 0.47, xi + 0.47, yi + 0.47, col = fill, border = NA)
  tc <- .on_color(fill)
  if (labels) graphics::text(xi, yi - if (values) 0.12 else 0, trt, col = tc, font = 2, cex = 0.85)
  if (values) graphics::text(xi, yi + if (labels) 0.2 else 0, formatC(v, format = "f", digits = digits), col = tc, cex = 0.75)
  graphics::axis(1, at = seq_along(cl), labels = cl, tick = FALSE)
  graphics::axis(2, at = seq_along(rl), labels = rl, tick = FALSE)
  .colorbar(pal, rng, digits)
  graphics::title(main = if (is.null(main)) paste0("Field map: ", if (value == "residuals") "residuals of " else "", fit$response) else main,
                  xlab = col, ylab = row, font.main = 1, cex.main = 1.05)
  invisible(data.frame(row = rr, col = cc, treatment = trt, value = v, stringsAsFactors = FALSE))
}

# ---------- residual diagnostics ----------

#' Residual plots of a design
#'
#' Six diagnostic plots of the residuals: residuals against fitted values,
#' normal Q-Q plot of the studentized residuals with a 95 % pointwise band,
#' scale-location, residuals by treatment, residuals in data order and a
#' histogram with the normal curve. Observations with |studentized residual|
#' above 2.5 are labelled with their row number.
#'
#' @param x An `"ad_anova"` object from [design_anova()].
#' @param which Plots to draw (1 to 6). Several plots share one page.
#' @param ... Not used.
#' @return Invisibly, a data frame with the fitted values, residuals and
#'   studentized residuals.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' plot(fit)
#' plot(fit, which = 2)
#' @export
plot.ad_anova <- function(x, which = 1:6, ...) {
  f <- x$qr_fit
  e <- f$residuals; fv <- x$y - e
  mse <- x$residual[["ms"]]
  stud <- e / (sqrt(mse) * sqrt(pmax(1e-12, 1 - .hat_values(f))))
  rown <- rownames(x$data)
  flag <- abs(stud) > 2.5
  grp <- .treatment_groups(x)
  col <- "#2f7d4f"
  mf <- if (length(which) == 1) c(1, 1) else if (length(which) == 2) c(1, 2) else if (length(which) <= 4) c(2, 2) else c(2, 3)
  op <- .style_par(mar = c(4.1, 4.4, 3, 1), mfrow = mf, oma = c(0, 0, 0.6, 0))
  on.exit(graphics::par(op))
  # row numbers of flagged points, to the left of points in the right half of the panel
  lab_at <- function(px, py, labs) {
    side <- ifelse(px > mean(graphics::par("usr")[1:2]), 2, 4)
    graphics::text(px, py, labs, pos = side, cex = 0.75, col = "#b5432f")
  }
  lab_out <- function(px, py) if (any(flag)) lab_at(px[flag], py[flag], rown[flag])
  for (w in which) {
    if (w == 1) {
      graphics::plot(fv, e, pch = 21, bg = col, col = "#ffffff", xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs fitted", font.main = 1)
      graphics::abline(h = 0, lty = 2, col = "#5f6368")
      if (length(unique(fv)) > 3) graphics::lines(stats::lowess(fv, e), col = "#b5432f", lwd = 1.5)
      lab_out(fv, e)
    } else if (w == 2) {
      n <- length(stud); o <- order(stud); th <- stats::qnorm(stats::ppoints(n))
      band <- cbind(stats::qnorm(stats::qbeta(0.025, seq_len(n), n + 1 - seq_len(n))),
                    stats::qnorm(stats::qbeta(0.975, seq_len(n), n + 1 - seq_len(n))))
      graphics::plot(th, stud[o], type = "n", xlab = "Normal quantiles", ylab = "Studentized residuals", main = "Normal Q-Q", font.main = 1,
                     ylim = range(c(stud, band)))
      graphics::polygon(c(th, rev(th)), c(band[, 1], rev(band[, 2])), col = "#2f7d4f22", border = NA)
      graphics::abline(0, 1, lty = 2, col = "#5f6368")
      graphics::points(th, stud[o], pch = 21, bg = col, col = "#ffffff")
      if (any(flag)) lab_at(th[match(which(flag), o)], stud[flag], rown[flag])
    } else if (w == 3) {
      graphics::plot(fv, sqrt(abs(stud)), pch = 21, bg = col, col = "#ffffff", xlab = "Fitted values",
                     ylab = expression(sqrt("|studentized residual|")), main = "Scale-location", font.main = 1)
      if (length(unique(fv)) > 3) graphics::lines(stats::lowess(fv, sqrt(abs(stud))), col = "#b5432f", lwd = 1.5)
    } else if (w == 4) {
      graphics::boxplot(e ~ grp, col = .lighten(col, 0.6), border = .darken(col, 0.3), xlab = "", ylab = "Residuals",
                        las = if (nlevels(factor(grp)) > 5) 2 else 1, cex.axis = 0.8)
      graphics::title(main = "Residuals by treatment", font.main = 1)
      graphics::abline(h = 0, lty = 2, col = "#5f6368")
    } else if (w == 5) {
      graphics::plot(seq_along(e), e, type = "b", pch = 21, bg = col, col = "#7f8c86", xlab = "Order in the data",
                     ylab = "Residuals", main = sprintf("Residuals in data order (DW = %.2f)", .durbin_watson(e)), font.main = 1)
      graphics::abline(h = 0, lty = 2, col = "#5f6368")
      lab_out(seq_along(e), e)
    } else if (w == 6) {
      h <- graphics::hist(e, plot = FALSE, breaks = "Sturges")
      s <- sqrt(sum(e^2) / x$residual[["df"]])
      xx <- seq(min(h$breaks), max(h$breaks), length.out = 200)
      graphics::hist(e, breaks = h$breaks, freq = FALSE, col = .lighten(col, 0.55), border = "#ffffff",
                     ylim = c(0, max(h$density, stats::dnorm(0, 0, s)) * 1.1), xlab = "Residuals", main = "Histogram of residuals", font.main = 1)
      graphics::lines(xx, stats::dnorm(xx, 0, s), col = "#b5432f", lwd = 2)
    }
  }
  invisible(data.frame(row = rown, fitted = fv, residual = e, studentized = stud))
}

#' @export
plot.ad_transformations <- function(x, main = "Box-Cox profile log-likelihood", ...) {
  b <- x$box_cox
  if (is.null(b)) stop("Box-Cox needs strictly positive values", call. = FALSE)
  op <- .style_par(); on.exit(graphics::par(op))
  graphics::plot(b$profile$lambda, b$profile$loglik, type = "l", lwd = 2.2, col = "#2f7d4f",
                 xlab = expression(lambda), ylab = "Log-likelihood", main = main, font.main = 1, ...)
  graphics::abline(h = b$cut, lty = 3, col = "#5f6368")
  graphics::abline(v = c(b$lower, b$upper), lty = 2, col = "#5f6368")
  graphics::abline(v = 1, lty = 3, col = "#b5432f")
  graphics::points(b$lambda, max(b$profile$loglik), pch = 21, bg = "#2f7d4f", col = "#ffffff", cex = 1.6)
  .note(sprintf("lambda = %.2f, 95 %% interval %.1f to %.1f", b$lambda, b$lower, b$upper))
  invisible(b$profile)
}

# ---------- export ----------

#' Save a figure to a file
#'
#' Opens the graphics device that matches the extension of `file` (`.png`,
#' `.tiff`, `.jpg`, `.pdf`, `.svg` or `.eps`), evaluates the plotting code
#' and closes the device. Raster formats carry the resolution in their
#' header, as journals ask.
#'
#' @param file Output file; its extension chooses the format.
#' @param expr Plotting code, for instance `plot(cmp)`; it is evaluated
#'   with the device open.
#' @param width,height Size of the figure.
#' @param units `"cm"`, `"mm"` or `"in"`.
#' @param res Resolution in dots per inch of raster formats.
#' @param pointsize Base size of the text, in points.
#' @return The path of the file, invisibly.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' f <- tempfile(fileext = ".png")
#' save_figure(f, plot(compare_means(fit, "Nitrogen_kg_ha")), width = 12, height = 9)
#' file.exists(f)
#' @export
save_figure <- function(file, expr, width = 16, height = 11, units = c("cm", "mm", "in"), res = 600, pointsize = 10) {
  units <- match.arg(units)
  k <- switch(units, cm = 1 / 2.54, mm = 1 / 25.4, `in` = 1)
  wi <- width * k; hi <- height * k
  ext <- tolower(sub(".*\\.", "", basename(file)))
  cairo <- isTRUE(capabilities("cairo"))
  switch(ext,
    png = grDevices::png(file, width = wi, height = hi, units = "in", res = res, pointsize = pointsize,
                         type = if (cairo) "cairo" else getOption("bitmapType")),
    tif = , tiff = grDevices::tiff(file, width = wi, height = hi, units = "in", res = res, pointsize = pointsize,
                                   compression = "lzw", type = if (cairo) "cairo" else getOption("bitmapType")),
    jpg = , jpeg = grDevices::jpeg(file, width = wi, height = hi, units = "in", res = res, pointsize = pointsize,
                                   quality = 95, type = if (cairo) "cairo" else getOption("bitmapType")),
    pdf = grDevices::pdf(file, width = wi, height = hi, pointsize = pointsize, useDingbats = FALSE),
    svg = { if (!cairo) stop("SVG output needs cairo support in this R", call. = FALSE)
            grDevices::svg(file, width = wi, height = hi, pointsize = pointsize) },
    eps = grDevices::postscript(file, width = wi, height = hi, pointsize = pointsize, horizontal = FALSE,
                                onefile = FALSE, paper = "special"),
    stop("unknown format '.", ext, "'; use png, tiff, jpg, pdf, svg or eps", call. = FALSE))
  dev <- grDevices::dev.cur()
  on.exit(if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev))
  force(expr)
  invisible(file)
}
