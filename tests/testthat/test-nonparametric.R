test_that("Friedman reproduces friedman.test and the figures of the web application", {
  f <- design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block")
  np <- nonparametric(f)
  expect_equal(np$method, "friedman")
  ft <- stats::friedman.test(Severity_score ~ Fungicide | Block, disease_scores)
  expect_equal(unname(np$test$statistic), unname(ft$statistic))
  expect_equal(round(unname(np$test$effect), 3), 0.936)
  expect_setequal(np$groups$group, c("a", "b", "c", "d"))
  expect_equal(np$groups$level[1], "Control")
  nm <- nonparametric(f, posthoc = "nemenyi")
  expect_equal(nm$posthoc, "Nemenyi")
  expect_equal(nm$p_adjust, "none")
})

test_that("Kruskal-Wallis, Welch and Mann-Whitney agree with base R", {
  f <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
  np <- nonparametric(f)
  expect_equal(np$method, "kruskal")
  kt <- stats::kruskal.test(Yield_kg_plot ~ Variety, bean_varieties)
  expect_equal(unname(np$test$statistic), unname(kt$statistic))
  expect_equal(round(np$test$p, 4), 0.0015)
  expect_equal(np$pairs$p_adj, stats::p.adjust(np$pairs$p, "holm"))
  w <- nonparametric(f, method = "welch")
  ow <- stats::oneway.test(Yield_kg_plot ~ Variety, bean_varieties)
  expect_equal(unname(w$test$statistic), unname(ow$statistic))
  mw <- nonparametric(f, method = "kruskal_mw")
  a <- bean_varieties$Yield_kg_plot[bean_varieties$Variety == mw$pairs$level_i[1]]
  b <- bean_varieties$Yield_kg_plot[bean_varieties$Variety == mw$pairs$level_j[1]]
  ref <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE, correct = FALSE))$p.value
  expect_equal(mw$pairs$p[1], ref)
})

test_that("methods that do not fit the design are refused", {
  f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "split_rcbd")
  expect_error(nonparametric(f), "no rank-based method")
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  expect_error(nonparametric(f, method = "kruskal"), "does not fit")
})

test_that("ART aligned sums are zero and SRH uses ranks", {
  set.seed(11)
  fd <- expand.grid(rep = 1:4, A = c("a1", "a2", "a3"), B = c("b1", "b2"))
  fd$y <- stats::rexp(nrow(fd), 1 / (3 + as.integer(fd$A)))
  f <- design_anova(fd, "y", c("A", "B"))
  art <- nonparametric(f, method = "art")
  expect_true(all(abs(art$table$aligned_sum) < 1e-8))
  srh <- nonparametric(f)
  expect_equal(srh$method, "srh")
  rk <- rank(fd$y)
  expect_equal(sum(srh$table$ss), sum((rk - mean(rk))^2) -
                 sum(stats::residuals(stats::lm(rk ~ A * B, fd))^2))
})

test_that("compact letters handle non-transitive results", {
  dm <- matrix(FALSE, 3, 3); dm[1, 3] <- dm[3, 1] <- TRUE
  expect_equal(AgriDesign:::.letters(c(3, 2, 1), dm), c("a", "ab", "b"))
  expect_equal(AgriDesign:::.letters(c(1, 2, 3), matrix(FALSE, 3, 3)), c("a", "a", "a"))
})
