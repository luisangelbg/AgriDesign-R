concurrence <- function(trt, blk) { N <- unclass(table(trt, blk)); M <- N %*% t(N); M[upper.tri(M)] }

test_that("Graeco-Latin squares are Latin in both letters and orthogonal", {
  for (t in c(3, 4, 5, 8)) {
    b <- make_design("graeco", trt = t, seed = t)$book
    expect_true(all(table(b$Row, b$Treatment) == 1) && all(table(b$Column, b$Greek) == 1))
    expect_true(all(table(b$Treatment, b$Greek) == 1))
  }
  expect_error(make_design("graeco", trt = 6, seed = 1), "none exists of order 2 or 6")
})

test_that("stored difference sets are difference sets", {
  for (key in names(AgriDesign:::.difference_sets)) {
    v <- as.integer(sub("-.*", "", key)); D <- AgriDesign:::.difference_sets[[key]]; k <- length(D)
    dif <- outer(D, D, "-") %% v; dif <- dif[row(dif) != col(dif)]
    expect_true(all(tabulate(dif, v - 1) == k * (k - 1) / (v - 1)), info = key)
  }
})

test_that("Youden squares and BIBDs are balanced", {
  b <- make_design("youden", trt = 7, k = 3, seed = 2)$book
  expect_true(all(table(b$Row, b$Treatment) == 1))
  expect_equal(unique(concurrence(b$Treatment, b$Column)), 1)
  b <- make_design("bib", trt = 6, k = 3, seed = 2)$book
  expect_equal(length(unique(concurrence(b$Treatment, b$Block))), 1)
  expect_equal(design_efficiency(b, "Treatment", "Block")$E, 4 * 6 / (10 * 3))
  expect_error(make_design("youden", trt = 9, k = 4, seed = 1), "no cyclic difference set")
})

test_that("efficiency factors match the closed forms", {
  lay <- make_design("lattice", trt = 16, k = 4, r = 2, seed = 1)
  expect_equal(design_efficiency(lay)$E, 5 / 7)
  expect_equal(design_efficiency(make_design("rcbd", trt = 5, r = 3, seed = 1))$E, 1)
  lay <- make_design("alpha", trt = 12, k = 3, r = 2, seed = 1)
  e <- design_efficiency(lay)
  expect_equal(e$upper_bound, 11 / (11 + 2 * 3))
  expect_lte(e$E, e$upper_bound + 1e-12)
  expect_equal(round(e$E, 4), 0.5946)
  expect_true(all(table(lay$book$Block, lay$book$Treatment) == 1))
})

test_that("cyclic designs are equireplicated and connected", {
  lay <- make_design("cyclic", trt = 10, k = 3, r = 6, seed = 1)
  expect_true(all(table(lay$book$Treatment) == 6))
  expect_true(design_efficiency(lay)$connected)
  expect_error(make_design("cyclic", trt = 10, k = 3, r = 4, seed = 1), "multiple of the block size")
})

test_that("new designs are suggested and analysed", {
  lay <- make_design("bib", trt = 7, k = 3, seed = 3)
  b <- field_book(lay, "y"); b$y <- 10 + as.integer(factor(b$Treatment)) / 3 + seq_len(nrow(b)) %% 3 / 5
  f <- design_anova(b, "y", "Treatment", block = "Block")
  expect_equal(f$design, "bib")
  s <- lay$skeleton
  expect_equal(c(f$table$df, f$residual[["df"]]), s$df[s$source != "Total"])
  lay <- make_design("graeco", trt = 5, seed = 3)
  b <- field_book(lay, "y"); b$y <- 5 + seq_len(nrow(b)) %% 4
  f <- design_anova(b, "y", "Treatment", block = "Greek", row = "Row", col = "Column")
  expect_equal(f$design, "graeco")
  expect_equal(f$residual[["df"]], (5 - 1) * (5 - 3))
  lay <- make_design("youden", trt = 7, k = 4, seed = 3)
  b <- field_book(lay, "y"); b$y <- 5 + seq_len(nrow(b)) %% 5
  expect_equal(design_anova(b, "y", "Treatment", row = "Row", col = "Column")$design, "youden")
})

test_that("REML recovery gives the GLS means of a random-block model", {
  lay <- make_design("alpha", trt = 12, k = 3, r = 3, seed = 4)
  b <- field_book(lay, "y")
  blk <- as.integer(factor(paste(b$Block, b$IncBlock)))
  b$y <- 10 + as.integer(factor(b$Treatment)) / 5 + c(-1.2, 0.4, 0.9, -0.3, 1.1, -0.8, 0.2, -0.5, 0.7, -1.0, 0.6, 0.3)[blk] +
    c(0.2, -0.1, 0.05, -0.2, 0.15)[seq_len(nrow(b)) %% 5 + 1]
  f <- design_anova(b, "y", "Treatment", block = c("Block", "IncBlock"))
  rc <- recover_interblock(f)
  expect_true(rc$sigma2_block > 0)
  expect_equal(nrow(rc$means), 12)
  expect_lte(rc$sed[["combined"]], rc$sed[["intra_block"]] + 1e-12)
  expect_error(recover_interblock(design_anova(bean_varieties, "Yield_kg_plot", "Variety")), "no blocks")
})
