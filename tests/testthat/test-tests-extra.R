test_that("Waller-Duncan critical t matches SAS and its limits", {
  # SAS/STAT ANOVA documentation (red clover): K = 100, q = 5, f = 24, F = 14.37
  expect_equal(round(waller_critical(100, 5, 24, 14.37), 5), 1.91873)
  # large F: the critical t no longer depends on q
  big <- waller_critical(100, c(1, 4, 10), 20, 1e8)
  expect_equal(big[1], big[3], tolerance = 1e-6)
  # f, F -> infinity: normal rule g(z) / g(-z) = K with g(z) = z Phi(z) + phi(z)
  g <- function(z) z * stats::pnorm(z) + stats::dnorm(z)
  z <- stats::uniroot(function(z) g(z) / g(-z) - 100, c(0.5, 5), tol = 1e-12)$root
  expect_equal(waller_critical(100, 3, 1e7, 1e8), z, tolerance = 1e-4)
  # smaller F, larger critical t; larger K, larger critical t
  tt <- waller_critical(100, 4, 20, c(2, 5, 50))
  expect_true(all(diff(tt) < 0))
  expect_gt(waller_critical(500, 4, 20, 5), waller_critical(100, 4, 20, 5))
  expect_error(waller_critical(1, 4, 20, 5), "greater than 1")
})

test_that("waller_duncan() returns a comparison that prints and plots", {
  fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  wd <- waller_duncan(fit, "Nitrogen_kg_ha")
  expect_s3_class(wd, "ad_comparison")
  expect_equal(wd$method, "waller")
  ref <- fit$table[fit$table$term == "Nitrogen_kg_ha", ]
  expect_equal(wd$t_critical, waller_critical(100, ref$df, ref$df_error, ref$F))
  se <- sqrt(2 * fit$residual[["ms"]] / 4)
  expect_equal(unname(wd$msd), wd$t_critical * se, tolerance = 1e-8)
  expect_equal(wd$pairs$different, abs(wd$pairs$diff) > wd$t_critical * se)
  expect_setequal(wd$means$group, c("a", "b", "c", "d"))
  expect_output(print(wd), "Waller-Duncan k-ratio t test \\(K = 100\\)")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(nrow(plot(wd)), 4)
  expect_equal(waller_duncan(fit, "Nitrogen_kg_ha", alpha = 0.01)$K, 500)
  expect_error(waller_duncan(fit, "Nitrogen_kg_ha", alpha = 0.2), "0.10, 0.05 or 0.01")
  expect_error(waller_duncan(fit, "Block"), "not a treatment factor")
})

test_that("Durbin's test: Friedman with complete blocks and balanced incomplete blocks", {
  f <- design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block")
  d <- durbin_test(f)
  ft <- stats::friedman.test(Severity_score ~ Fungicide | Block, disease_scores)
  expect_equal(unname(d$test$statistic), unname(ft$statistic))
  du <- durbin_test(bib_sorghum, "Yield_t_ha", "Hybrid", "Block", p_adjust = "none")
  expect_s3_class(du, c("ad_durbin", "ad_nonparametric"))
  expect_equal(unname(du$design[c("treatments", "blocks", "k", "r", "lambda")]), c(7, 7, 3, 3, 1))
  # hand computation of T1 = (t - 1) sum (R_j - r (k + 1) / 2)^2 / (A - C)
  rk <- stats::ave(bib_sorghum$Yield_t_ha, bib_sorghum$Block, FUN = rank)
  R <- tapply(rk, bib_sorghum$Hybrid, sum)
  T1 <- 6 * sum((R - 3 * 2)^2) / (sum(rk^2) - 7 * 3 * 16 / 4)
  expect_equal(unname(du$test$statistic), unname(T1))
  expect_equal(sum(du$groups$rank_sum), 7 * 6)
  expect_output(AgriDesign:::print.ad_durbin(du), "Least significant difference of rank sums")
  bad <- bib_sorghum[-1, ]
  expect_error(durbin_test(bad, "Yield_t_ha", "Hybrid", "Block"), "equal size")
  expect_error(durbin_test(bean_varieties, "Yield_kg_plot", "Variety"), "needs a block")
})

test_that("van der Waerden scores and pairwise tests", {
  vw <- waerden_test(bean_varieties, "Yield_kg_plot", "Variety", p_adjust = "none")
  y <- bean_varieties$Yield_kg_plot; g <- factor(bean_varieties$Variety, levels = unique(bean_varieties$Variety))
  A <- stats::qnorm(rank(y) / (length(y) + 1))
  T1 <- sum(table(g) * tapply(A, g, mean)^2) / (sum(A^2) / (length(y) - 1))
  expect_equal(unname(vw$test$statistic), unname(T1))
  expect_equal(sort(vw$groups$mean_score), sort(as.numeric(tapply(A, g, mean))))
  expect_equal(nrow(vw$pairs), 10)
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  expect_error(waerden_test(f), "without blocks")
  expect_error(waerden_test(bean_varieties, "Yield_kg_plot"), "give 'response' and 'treatment'")
})

test_that("Mood's median test agrees with chisq.test and fisher.test", {
  y <- bean_varieties$Yield_kg_plot; g <- factor(bean_varieties$Variety, levels = unique(bean_varieties$Variety))
  mt <- median_test(bean_varieties, "Yield_kg_plot", "Variety", p_adjust = "none")
  ref <- suppressWarnings(stats::chisq.test(table(y > stats::median(y), g), correct = FALSE))
  expect_equal(unname(mt$test$statistic), unname(ref$statistic))
  expect_equal(mt$test$p, ref$p.value)
  s <- g %in% levels(g)[1:2]; yy <- y[s]
  p2 <- suppressWarnings(stats::chisq.test(table(yy > stats::median(yy), droplevels(g[s])), correct = TRUE)$p.value)
  expect_equal(median_test(bean_varieties, "Yield_kg_plot", "Variety", correct = TRUE, p_adjust = "none")$pairs$p[1], p2)
  pe <- stats::fisher.test(table(yy > stats::median(yy), droplevels(g[s])))$p.value
  expect_equal(median_test(bean_varieties, "Yield_kg_plot", "Variety", exact = TRUE)$pairs$p[1], pe)
  expect_equal(mt$table$above + mt$table$at_or_below, mt$table$n)
  expect_output(AgriDesign:::print.ad_median_test(mt), "Grand median")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(nrow(plot(mt)), 5)
})
