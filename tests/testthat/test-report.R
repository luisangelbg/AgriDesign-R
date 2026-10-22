test_that("the report prints the numbers of the analysis", {
  fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  h <- design_report(fit, figures = FALSE, units = "t/ha")
  expect_true(grepl(">437.267<", h, fixed = TRUE))
  expect_true(grepl("Minimum significant difference = 0.283", h, fixed = TRUE))
  expect_true(grepl("randomised complete block design", h, fixed = TRUE))
  expect_true(grepl("<span class=\"grade\">A</span>", h, fixed = TRUE))
  expect_false(grepl("&amp;[a-z]+;", h))
  expect_false(grepl("<script|<link", h))
  expect_equal(lengths(gregexpr("<table>", h, fixed = TRUE)), attr(h, "tables"))
})

test_that("sections can be left out and the file is written", {
  fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
  h <- design_report(fit, figures = FALSE, assumptions = FALSE, pairwise = FALSE, cite = FALSE)
  expect_false(grepl("Assumptions of the analysis", h, fixed = TRUE))
  expect_false(grepl("Appendix A", h, fixed = TRUE))
  expect_false(grepl("How to cite", h, fixed = TRUE))
  f <- tempfile(fileext = ".html")
  design_report(fit, f, figures = FALSE, data = TRUE, method = "dunnett", control = "Negro Jamapa")
  h <- paste(readLines(f, encoding = "UTF-8"), collapse = "\n")
  expect_true(grepl("Appendix B", h, fixed = TRUE))
  expect_true(grepl("differs from the control", h, fixed = TRUE))
})

test_that("figures are inlined with their own ids", {
  fit <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd")
  h <- design_report(fit, assumptions = FALSE, pairwise = FALSE)
  expect_gte(attr(h, "figures"), 3)
  ids <- regmatches(h, gregexpr("id=\"[^\"]+\"", h))[[1]]
  expect_equal(anyDuplicated(ids), 0)
  expect_true(grepl("pooled error", h, fixed = TRUE))
})

test_that("the base64 encoder matches the standard alphabet", {
  expect_equal(AgriDesign:::.base64(charToRaw("Man")), "TWFu")
  expect_equal(AgriDesign:::.base64(charToRaw("Ma")), "TWE=")
  expect_equal(AgriDesign:::.base64(charToRaw("M")), "TQ==")
  expect_equal(AgriDesign:::.base64(charToRaw("AgriDesign")), "QWdyaURlc2lnbg==")
})
