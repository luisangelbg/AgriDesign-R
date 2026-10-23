# Draws the hexagon logo of AgriDesign with base R graphics (no images or fonts
# of others): a field trial in perspective, laid out as a Latin square, with a
# maize plant growing from it, on the forest green and gold of AgriDesign.
# Run from the package root:  Rscript data-raw/make-logo.R
# Writes man/figures/logo.png (pkgdown and README), logo.svg and a 600 dpi PNG.

gold <- "#c9a45c"; gold_light <- "#e6cf98"; green_dark <- "#16382a"; green <- "#1f4d33"
cream <- "#f6f1e3"
plot_cols <- c("#5aa469", "#c8842a", "#3f86c0", "#b5432f")   # the treatments of the square

# pointy-top hexagon of width w centred at (0, 0)
hex <- function(w) { a <- pi / 6 + (0:5) * pi / 3; r <- w / sqrt(3); cbind(r * cos(a + pi / 3), r * sin(a + pi / 3)) }

draw_logo <- function() {
  op <- par(mar = c(0, 0, 0, 0), bg = "transparent"); on.exit(par(op))
  plot.new(); plot.window(xlim = c(-1, 1), ylim = c(-1.1547, 1.1547), asp = 1, xaxs = "i", yaxs = "i")
  outer <- hex(2); inner <- hex(1.86)
  polygon(outer, col = gold, border = NA)
  polygon(inner, col = green, border = NA)
  # maize plant (drawn first, so that it grows from behind the field): stalk, leaves, ear, tassel
  base_y <- 0.02; top_y <- 0.80
  polygon(c(-0.024, 0.024, 0.018, -0.018), c(base_y, base_y, top_y, top_y), col = gold_light, border = NA)
  leaf <- function(y, dir, len, lift) {
    t <- seq(0, 1, length.out = 60)
    xs <- dir * (0.02 + len * t)
    yc <- y + lift * sin(pi * t * 0.8) - 0.16 * t^2       # rises from the stalk and droops at the tip
    w <- 0.05 * sin(pi * t^0.8) * (1 - 0.3 * t)
    polygon(c(xs, rev(xs)), c(yc + w, rev(yc - w)), col = gold_light, border = NA)
  }
  leaf(0.20, -1, 0.46, 0.22); leaf(0.30, 1, 0.44, 0.24); leaf(0.46, -1, 0.34, 0.18); leaf(0.56, 1, 0.30, 0.14)
  ex <- 0.075; ey <- 0.44; t <- seq(0, 2 * pi, length.out = 60)
  polygon(ex + 0.05 * cos(t), ey + 0.12 * sin(t), col = gold, border = NA)
  for (k in -2:2) points(ex + 0.022 * c(-1, 1), rep(ey + k * 0.042, 2), pch = 16, cex = 0.32 * .cexk, col = "#8f6d2c")
  for (a in c(-0.7, -0.35, 0, 0.35, 0.7)) segments(0, top_y, 0.13 * sin(a), top_y + 0.11 * cos(a), col = gold_light, lwd = 2.4 * .cexk)
  # field in perspective: a Latin square of 4 x 4 plots mapped onto a trapezoid
  y0 <- -0.40; y1 <- 0.06; x0b <- 0.68; x0t <- 0.40
  to_field <- function(u, v) cbind((-1 + 2 * u) * (x0b + (x0t - x0b) * v), y0 + (y1 - y0) * v)
  latin <- matrix(c(1, 2, 3, 4, 2, 4, 1, 3, 3, 1, 4, 2, 4, 3, 2, 1), 4, byrow = TRUE)
  g <- 0.04
  for (i in 1:4) for (j in 1:4) {
    u <- c(j - 1 + g, j - g, j - g, j - 1 + g) / 4
    v <- c(i - 1 + g, i - 1 + g, i - g, i - g) / 4
    polygon(to_field(u, v), col = plot_cols[latin[i, j]], border = NA)
  }
  # name, scaled to a fixed share of the hexagon whatever the resolution
  cx <- 1.32 / strwidth("AgriDesign", font = 2, cex = 1, family = "sans")
  text(0, -0.63, "AgriDesign", col = cream, font = 2, cex = cx, family = "sans")
}

# line widths and point sizes follow the size of the drawing (cexk = 1 at 278 px of height)
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
w_in <- 1.73; h_in <- 2
.cexk <- 1
png("man/figures/logo.png", width = 240, height = 278, bg = "transparent", type = "cairo", res = 139)
draw_logo(); dev.off()
.cexk <- 1
png("man/figures/logo-600dpi.png", width = w_in * 600, height = h_in * 600, bg = "transparent", type = "cairo", res = 600)
draw_logo(); dev.off()
svg("man/figures/logo.svg", width = w_in, height = h_in, bg = "transparent")
draw_logo(); dev.off()
