fac <- function(d, v) { for (x in v) d[[x]] <- factor(d[[x]]); d }
row_of <- function(fit, term, col) fit$table[fit$table$term == term, col]

test_that("RCBD reproduces the figures of the web application", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  expect_equal(f$design, "rcbd")
  expect_equal(row_of(f, "Nitrogen_kg_ha", "df"), 3)
  expect_equal(f$residual[["df"]], 9)
  expect_equal(round(row_of(f, "Nitrogen_kg_ha", "F"), 3), 437.267)
  expect_equal(round(f$cv, 1), 2.1)
  m <- ls_means(f, "Nitrogen_kg_ha")
  expect_equal(m$mean, c(4.215, 5.845, 6.885, 7.1825), tolerance = 1e-4)
  expect_equal(m$se, rep(sqrt(f$residual[["ms"]] / 4), 4))
})

test_that("CRD, Latin square and factorial agree with aov()", {
  f <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
  a <- summary(stats::aov(Yield_kg_plot ~ Variety, bean_varieties))[[1]]
  expect_equal(row_of(f, "Variety", "F"), a[1, "F value"])

  f <- design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column")
  expect_equal(f$design, "latin")
  a <- summary(stats::aov(Yield_t_ha ~ Row + Column + Treatment, fac(wheat_latin, c("Row", "Column"))))[[1]]
  expect_equal(row_of(f, "Treatment", "F"), a[3, "F value"])
  expect_equal(round(row_of(f, "Treatment", "F"), 2), 164.10)

  f <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
                    block = "Block", ss_type = 1)
  a <- summary(stats::aov(Yield_kg_plant ~ Block + Variety * Fertilizer, fac(tomato_factorial, "Block")))[[1]]
  expect_equal(f$table$ss, a[1:4, "Sum Sq"])
})

test_that("split-plot and strip-plot use the right error strata", {
  d <- fac(split_irrigation, "Block")
  f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "split_rcbd")
  expect_equal(row_of(f, "Irrigation", "error"), "Block:Irrigation")
  expect_equal(row_of(f, "Irrigation", "df_error"), 2)
  expect_equal(round(row_of(f, "Irrigation", "F"), 1), 10460.6)
  expect_true(is.na(row_of(f, "Block:Irrigation", "F")))

  f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "split_crd")
  expect_equal(row_of(f, "Irrigation:Block", "df"), 4)
  a <- suppressWarnings(summary(stats::aov(Yield_t_ha ~ Irrigation * Variety + Error(Irrigation:Block), d)))
  expect_equal(row_of(f, "Irrigation", "F"), a[["Error: Irrigation:Block"]][[1]]["Irrigation", "F value"])

  f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "strip_rcbd")
  a <- suppressWarnings(summary(stats::aov(
    Yield_t_ha ~ Block + Irrigation * Variety + Error(Block:Irrigation + Block:Variety), d)))
  expect_equal(row_of(f, "Variety", "F"), a[["Error: Block:Variety"]][[1]]["Variety", "F value"])
  expect_equal(row_of(f, "Irrigation:Variety", "F"),
               a[["Error: Within"]][[1]]["Irrigation:Variety", "F value"])
})

test_that("nested design codes the inner factor within the outer one", {
  set.seed(1)
  nd <- expand.grid(rep = 1:3, Plant = c("p1", "p2", "p3"), Orchard = c("O1", "O2", "O3", "O4"))
  nd$y <- stats::rnorm(nrow(nd), 10 + as.integer(nd$Orchard))
  f <- design_anova(nd, "y", c("Orchard", "Plant"), design = "nested")
  expect_equal(row_of(f, "Orchard:Plant", "df"), 8)
  a <- suppressWarnings(summary(stats::aov(y ~ Orchard + Error(Orchard:Plant), nd)))
  expect_equal(row_of(f, "Orchard", "F"), a[["Error: Orchard:Plant"]][[1]]["Orchard", "F value"])
})

test_that("ANCOVA adjusts the means to the mean of the covariate", {
  d <- maize_nitrogen
  d$Stand <- c(58, 61, 63, 57, 60, 62, 59, 64, 55, 61, 60, 62, 63, 58, 59, 60)
  f <- design_anova(d, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block", covariate = "Stand")
  dd <- fac(d, c("Block", "Nitrogen_kg_ha"))
  fit <- stats::lm(Yield_t_ha ~ Stand + Block + Nitrogen_kg_ha, dd)
  nd <- data.frame(Stand = mean(d$Stand), Block = factor(rep(levels(dd$Block), 4)),
                   Nitrogen_kg_ha = factor(rep(levels(dd$Nitrogen_kg_ha), each = 4),
                                           levels = levels(dd$Nitrogen_kg_ha)))
  ref <- tapply(stats::predict(fit, nd), nd$Nitrogen_kg_ha, mean)
  expect_equal(ls_means(f, "Nitrogen_kg_ha")$mean, unname(as.numeric(ref)))
})

test_that("roles that do not fit the design are refused with a clear message", {
  expect_error(design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", design = "rcbd"),
               "needs 1 block column")
  expect_error(design_anova(maize_nitrogen, "Yield", "Nitrogen_kg_ha"), "not found")
  expect_error(design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", design = "xyz"),
               "unknown design")
  expect_equal(nrow(designs()), 13)
})
