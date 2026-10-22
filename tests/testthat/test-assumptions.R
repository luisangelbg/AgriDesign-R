tst <- function(a, name, col = "p") a$tests[a$tests$test == name, col]

test_that("maize checks reproduce the figures of the web application", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  a <- check_assumptions(f)
  expect_equal(a$grade, "A")
  expect_equal(round(tst(a, "Shapiro-Wilk"), 4), 0.2477)
  expect_equal(round(tst(a, "Levene (median)"), 4), 0.6672)
  expect_equal(round(tst(a, "Levene (mean)"), 4), 0.0287)
  expect_equal(round(tst(a, "Tukey one-degree-of-freedom"), 4), 0.8208)
  expect_equal(round(tst(a, "Anderson-Darling"), 4), 0.1840)
  expect_equal(a$suspicious$row, 16)
  e <- a$residuals$residual
  expect_equal(tst(a, "Shapiro-Wilk", "statistic"), unname(stats::shapiro.test(e)$statistic))
  fit <- stats::lm(Yield_t_ha ~ factor(Block) + factor(Nitrogen_kg_ha), maize_nitrogen)
  expect_equal(a$residuals$studentized, unname(stats::rstandard(fit)))
  bc <- box_cox(f)
  expect_equal(c(bc$lambda, bc$lower, bc$upper), c(0.7, 0, 1.5))
  expect_true(bc$contains_one)
})

test_that("Tukey's additivity test equals the textbook formula for an RCBD", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  ym <- tapply(maize_nitrogen$Yield_t_ha, list(maize_nitrogen$Block, maize_nitrogen$Nitrogen_kg_ha), mean)
  bi <- rowMeans(ym) - mean(ym); tj <- colMeans(ym) - mean(ym)
  ssn <- sum(ym * outer(bi, tj))^2 / (sum(bi^2) * sum(tj^2))
  expect_equal(AgriDesign:::.tukey_additivity(f$X, f$y)[["ss"]], ssn)
})

test_that("tomato fails additivity and the logarithm fixes it", {
  f <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block")
  a <- check_assumptions(f)
  expect_equal(a$grade, "C")
  expect_equal(round(tst(a, "Tukey one-degree-of-freedom", "statistic"), 2), 54.25)
  tr <- transformations(f)
  expect_equal(tr$suggested, "log")
  expect_equal(c(tr$box_cox$lambda, tr$box_cox$lower, tr$box_cox$upper), c(0.2, 0.1, 0.4))
  d <- tomato_factorial; d$lnY <- log(d$Yield_kg_plant)
  expect_equal(check_assumptions(design_anova(d, "lnY", c("Variety", "Fertilizer"), block = "Block"))$grade, "A")
})

test_that("a CRD has no additivity test and no transformation is suggested when none is needed", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha")
  a <- check_assumptions(f)
  expect_false(a$has_additivity)
  tr <- transformations(design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block"))
  expect_null(tr$suggested)
  expect_match(tr$verdict, "No transformation needed")
})

test_that("Box-Cox refuses zeros and transformations follow the range of the data", {
  d <- maize_nitrogen; d$Yield_t_ha[1] <- 0
  expect_null(box_cox(design_anova(d, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")))
  ids <- names(AgriDesign:::.transform_list(c(0, 0.2, 0.5)))
  expect_true(all(c("asinp", "sqrt", "log1") %in% ids))
  expect_false("log" %in% ids)
})
