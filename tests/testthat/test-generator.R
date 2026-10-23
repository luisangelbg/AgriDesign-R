test_that("layouts have the structure of their design", {
  b <- make_design("rcbd", trt = 5, r = 4, seed = 11)$book
  expect_true(all(table(b$Block, b$Treatment) == 1))
  b <- make_design("latin", trt = LETTERS[1:6], seed = 11)$book
  expect_true(all(table(b$Row, b$Treatment) == 1) && all(table(b$Column, b$Treatment) == 1))
  b <- make_design("split_rcbd", A = c("a1", "a2", "a3"), B = c("b1", "b2"), r = 2, seed = 11)$book
  expect_true(all(tapply(b$FieldCol, list(b$Block, b$A), function(v) length(unique(v))) == 1))
  b <- make_design("strip", A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 2, seed = 11)$book
  expect_true(all(tapply(b$FieldRow, list(b$Block, b$A), function(v) length(unique(v))) == 1))
  b <- make_design("ibd", trt = 9, r = 2, k = 3, seed = 11)$book
  expect_true(all(table(paste(b$Block, b$IncBlock)) == 3))
  b <- make_design("augmented", checks = c("C1", "C2"), entries = paste0("E", 1:7), r = 3, seed = 11)$book
  expect_equal(sum(b$Type == "entry"), 7)
  expect_true(all(table(b$Block[b$Type == "check"]) == 2))
  b <- make_design("fact_crd", A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 2, seed = 11)$book
  expect_equal(nrow(b), 12)
  expect_true(all(b$Treatment == paste(b$A, b$B, sep = " x ")))
})

test_that("plot numbers follow the serpentine path and the block option", {
  b <- make_design("rcbd", trt = 4, r = 3, seed = 2)$book
  b <- b[order(b$Plot), ]
  expect_equal(b$Plot, 1:12)
  expect_true(all(abs(diff(b$FieldRow)) + abs(diff(b$FieldCol)) == 1))
  b <- make_design("rcbd", trt = 4, r = 3, seed = 2, numbering = "block")$book
  expect_true(all(b$Plot %/% 100 == b$Block))
})

test_that("the seed reproduces the layout and leaves the session random state alone", {
  expect_identical(make_design("rcbd", trt = 6, r = 3, seed = 5)$book, make_design("rcbd", trt = 6, r = 3, seed = 5)$book)
  set.seed(42); a <- stats::runif(2)
  set.seed(42); invisible(make_design("latin", trt = 5, seed = 9)); b <- stats::runif(2)
  expect_identical(a, b)
  expect_true(is.numeric(make_design("crd", trt = 3, r = 2)$seed))
})

test_that("Jacobson-Matthews squares are Latin and reach every square of order 4", {
  is_latin <- function(m) all(apply(m, 1, anyDuplicated) == 0) && all(apply(m, 2, anyDuplicated) == 0)
  for (t in c(2, 3, 5, 7, 10)) expect_true(is_latin(AgriDesign:::.with_seed(t, AgriDesign:::.latin_square(t))))
  # Permuting a cyclic square of order 4 gives only the 432 squares with 4
  # intercalates (2 x 2 Latin subsquares); the other 144 have 12. A uniform
  # draw gives these with probability 1/4 (the full uniformity test is in
  # validation/v05-generator.R).
  intercalates <- function(m) {
    n <- 0
    for (i in 1:3) for (k in (i + 1):4) for (j in 1:3) for (l in (j + 1):4)
      n <- n + (m[i, j] == m[k, l] && m[i, l] == m[k, j])
    n
  }
  ni <- vapply(1:80, function(s) intercalates(AgriDesign:::.with_seed(s, AgriDesign:::.latin_square(4))), 0)
  expect_setequal(unique(ni), c(4, 12))
  expect_gt(mean(ni == 12), 0.1)
})

test_that("the field book goes straight into design_anova with the skeleton df", {
  lay <- make_design("split_rcbd", A = c("I1", "I2"), B = paste0("V", 1:4), r = 3, seed = 1)
  b <- field_book(lay, "Yield")
  expect_true(all(is.na(b$Yield)))
  set.seed(3); b$Yield <- stats::rnorm(nrow(b), 10)
  f <- design_anova(b, "Yield", c("A", "B"), block = "Block", design = "split_rcbd")
  s <- lay$skeleton
  expect_equal(c(f$table$df, f$residual[["df"]]), s$df[s$source != "Total"])
  f <- tempfile(fileext = ".csv")
  field_book(lay, c("Yield", "Height"), file = f)
  expect_equal(names(utils::read.csv(f))[c(1, ncol(utils::read.csv(f)))], c("Plot", "Height"))
})

test_that("power uses the noncentral t and agrees with power.t.test", {
  p <- power_design(cv = 12, diff = 15, t = 2, design = "crd", r = 5)
  expect_equal(p$table$power, stats::power.t.test(n = 5, delta = 15, sd = 12, strict = TRUE)$power)
  expect_equal(power_design(cv = 15, diff = 15, t = 4)$replicates, 17)
  expect_equal(power_design(cv = 10, diff = 15, t = 4)$replicates, 8)
  expect_true(is.na(power_design(cv = 40, diff = 2, t = 4)$replicates))
})

test_that("bad requests are refused with clear messages", {
  expect_error(make_design("ibd", trt = 10, r = 2, k = 3), "multiple of the block size")
  expect_error(make_design("split_rcbd", A = "a1", B = c("b1", "b2")), "levels of 'A'")
  expect_error(make_design("augmented", checks = c("C1", "C2"), entries = c("C1", "E2")), "different names")
  expect_error(make_design("rcbd", trt = 4, r = 1), "at least 2")
})

test_that("designs('layout') lists every generated design and its analysis runs on the field book", {
  lay_tab <- designs("layout")
  expect_equal(nrow(lay_tab), 16)
  expect_true(all(lay_tab$analyse_as %in% designs()$id))
  args <- list(crd = list(trt = 4), rcbd = list(trt = 4), latin = list(trt = 4),
               fact_crd = list(A = c("a1", "a2"), B = c("b1", "b2", "b3")),
               fact_rcbd = list(A = c("a1", "a2"), B = c("b1", "b2"), C = c("c1", "c2")),
               split_rcbd = list(A = c("a1", "a2"), B = c("b1", "b2", "b3")),
               strip = list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3")),
               splitsplit = list(A = c("a1", "a2"), B = c("b1", "b2"), C = c("c1", "c2")),
               ibd = list(trt = 6, k = 3), augmented = list(checks = c("C1", "C2", "C3"), entries = paste0("E", 1:8)),
               graeco = list(trt = 5), youden = list(trt = 7, k = 3), bib = list(trt = 7, k = 3),
               cyclic = list(trt = 6, k = 3, r = 3), alpha = list(trt = 12, k = 3), lattice = list(trt = 9, k = 3))
  expect_setequal(names(args), lay_tab$id)
  for (d in names(args)) {
    lay <- do.call(make_design, c(list(d, seed = 11), args[[d]]))
    data <- lay$book
    data$response <- stats::rnorm(nrow(data), 10)
    fit <- eval(parse(text = lay$analysis))
    expect_s3_class(fit, "ad_anova")
    expect_equal(fit$design, lay_tab$analyse_as[lay_tab$id == d], info = d)
  }
  expect_output(print(make_design("rcbd", trt = 3, r = 2, seed = 1)), 'design = "rcbd"')
})
