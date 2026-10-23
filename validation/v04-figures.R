# Stage 4: renders every figure to PNG (validation/figures/) for visual review
# and checks that what is drawn equals the numbers of the analysis.
suppressMessages(devtools::load_all(quiet = TRUE))
out <- "validation/figures"; dir.create(out, showWarnings = FALSE)
ok <- 0; bad <- 0
chk <- function(what, got, ref, tol = 1e-9) {
  got <- as.numeric(got); ref <- as.numeric(ref)
  pass <- length(got) == length(ref) && max(abs(got - ref)) < tol
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %s\n", if (pass) "OK" else "FAIL", what))
}
fig <- function(name, expr, w = 16, h = 11) save_figure(file.path(out, paste0(name, ".png")), expr, width = w, height = h, res = 150)

mz <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
cmp <- compare_means(mz, "Nitrogen_kg_ha")
fig("01-bars", r <- plot(cmp))
chk("bars: heights = LS means", r$mean, cmp$means$mean)
chk("bars: error bars = mean +/- SE", r$upper - r$mean, cmp$means$se)
fig("02-points-ci", r <- plot(cmp, style = "points", error = "ci", values = TRUE))
chk("points: CI half-width = t(0.975, df) SE", r$upper - r$mean, qt(0.975, mz$residual[["df"]]) * cmp$means$se)
fig("03-box", plot(cmp, style = "box"))
fig("04-strip", plot(cmp, style = "strip"))
fig("05-residuals", r <- plot(mz), w = 24, h = 15)
chk("residual plots: studentized = rstandard", r$studentized,
    rstandard(lm(Yield_t_ha ~ factor(Block) + factor(Nitrogen_kg_ha), maize_nitrogen)))
tom <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block")
fig("06-boxcox", plot(transformations(tom)))
fig("07-trend", plot(trend(mz, "Nitrogen_kg_ha")))
fig("08-interaction", r <- plot_interaction(tom, "Fertilizer", "Variety"))
chk("interaction: cell means = ls_means", r$mean, ls_means(tom, c("Fertilizer", "Variety"))$mean)
fig("09-cells", plot_cell_means(tom, "Variety", "Fertilizer"))
sp <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd")
se <- compare_means(sp, "Variety", by = "Irrigation")
fig("10-simple-effects", r <- plot(se))
chk("grouped bars: heights = cell means", r$mean, se$means$mean)
fig("11-ss-partition", r <- plot_ss_partition(sp), h = 9)
chk("SS partition: percentages add to 100", sum(r$pct), 100)
chk("SS partition: total = total SS", sum(r$ss), sp$total[["ss"]], 1e-9)
bn <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
fig("12-differences", r <- plot_differences(compare_means(bn, "Variety")))
fig("13-dunnett", plot(compare_means(bn, "Variety", method = "dunnett", control = "Negro Jamapa")))
fig("14-nonparametric", plot(nonparametric(design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block"))))
lt <- design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column")
fig("15-field", r <- plot_field(lt), w = 14)
chk("field map: one cell per plot", nrow(unique(r[c("row", "col")])), 16)
fig("16-field-residuals", r <- plot_field(lt, value = "residuals"), w = 14)
chk("field map: residuals of the model", r$value, lt$qr_fit$residuals)
# every format
for (e in c("tiff", "jpg", "pdf", "svg", "eps")) {
  f <- file.path(tempdir(), paste0("t.", e))
  ok_e <- tryCatch({ save_figure(f, plot(cmp), width = 8, height = 6); file.exists(f) && file.size(f) > 1000 }, error = function(err) FALSE)
  chk(paste("save_figure writes", e), ok_e, TRUE)
}
# resolution in the PNG header (pHYs chunk: pixels per metre)
f <- file.path(tempdir(), "dpi.png"); save_figure(f, plot(cmp), width = 8, height = 6, res = 600)
b <- readBin(f, "raw", file.size(f)); i <- grepRaw("pHYs", b)
ppm <- sum(as.integer(b[(i + 4):(i + 7)]) * 256^(3:0))
chk("PNG carries 600 dpi in its header", round(ppm * 0.0254), 600)
chk("no device left open", length(dev.list()), 0)
cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
