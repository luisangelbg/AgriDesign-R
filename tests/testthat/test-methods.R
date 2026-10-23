test_that("residuals, fitted, nobs and anova behave as for a linear model", {
  fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  d <- maize_nitrogen; d$Block <- factor(d$Block); d$N <- factor(d$Nitrogen_kg_ha)
  lmf <- stats::lm(Yield_t_ha ~ Block + N, d)
  expect_equal(unname(residuals(fit)), unname(stats::residuals(lmf)))
  expect_equal(unname(fitted(fit)), unname(stats::fitted(lmf)))
  expect_equal(unname(residuals(fit, type = "studentized")), unname(stats::rstandard(lmf)))
  expect_equal(nobs(fit), 16)
  a <- anova(fit)
  expect_s3_class(a, "anova")
  expect_equal(a["Nitrogen_kg_ha", "F value"], fit$table$F[2])
  expect_equal(a["Residuals", "Df"], 9)
})

test_that("anova() reports the error stratum of each term", {
  sp <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd")
  a <- anova(sp)
  expect_equal(a["Irrigation", "Tested against"], "Block:Irrigation")
  expect_equal(a["Variety", "Tested against"], "Residuals")
})

test_that("as.data.frame returns the means with their letters", {
  fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  m <- as.data.frame(compare_means(fit, "Nitrogen_kg_ha"))
  expect_true(all(c("mean", "se", "group") %in% names(m)))
  expect_equal(nrow(m), 4)
})
