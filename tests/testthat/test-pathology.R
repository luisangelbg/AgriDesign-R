test_that("AUDPC and AUDPS follow their formulas", {
  t <- c(0, 7, 14, 28); y <- c(2, 10, 30, 60)
  a <- sum((y[-1] + y[-4]) / 2 * diff(t))
  expect_equal(audpc(y, t), a)
  expect_equal(audps(y, t), a + (2 + 60) / 2 * 28 / 3)
  expect_equal(audpc(y, t, "relative"), a / (100 * 28))
  expect_equal(audps(y, t, "relative"), audps(y, t) / (100 * 28 * 4 / 3))
  expect_equal(audps(y, t, "standardized"), audps(y, t) / (28 * 4 / 3))
  # equal intervals: AUDPS = sum(y) * interval; the maximum relative areas are 1
  expect_equal(audps(c(5, 9, 20), c(10, 20, 30)), 34 * 10)
  expect_equal(audps(rep(100, 6), 1:6, "relative"), 1)
  expect_equal(audpc(rep(1, 3), c(1, 4, 9), "relative", ymax = 1), 1)
  # several units, missing assessments, Date values
  m <- rbind(p1 = y, p2 = c(1, NA, 20, 40))
  expect_equal(names(audpc(m, t)), c("p1", "p2"))
  expect_equal(audpc(m, t)[["p2"]], audpc(c(1, 20, 40), t[-2]))
  expect_equal(audpc(y, as.Date("2026-06-01") + t), a)
  expect_error(audpc(y, c(0, 7, 7, 28)), "strictly increasing")
  expect_error(audpc(y, t[-1]), "one date per assessment")
  expect_error(audpc(-y, t), "negative")
  expect_warning(audpc(y * 2, t), "exceed")
})

test_that("audpc_by_unit gives one row per plot ready for design_anova()", {
  a <- audpc_by_unit(blight_progress, "Severity_pct", "Day", c("Block", "Treatment"))
  expect_equal(nrow(a), 16)
  expect_true(all(a$n == 7))
  one <- blight_progress[blight_progress$Block == 2 & blight_progress$Treatment == "Systemic", ]
  i <- which(a$Block == 2 & a$Treatment == "Systemic")
  expect_equal(a$audps[i], audps(one$Severity_pct, one$Day))
  shuffled <- blight_progress[rev(seq_len(nrow(blight_progress))), ]
  b <- audpc_by_unit(shuffled, "Severity_pct", "Day", c("Block", "Treatment"))
  expect_equal(b$audpc[match(paste(a$Block, a$Treatment), paste(b$Block, b$Treatment))], a$audpc)
  fit <- design_anova(a, "audpc_rel", "Treatment", block = "Block")
  expect_lt(fit$table$p[fit$table$term == "Treatment"], 0.001)
  dup <- rbind(blight_progress, blight_progress[1, ])
  expect_error(audpc_by_unit(dup, "Severity_pct", "Day", c("Block", "Treatment")), "same time")
  short <- blight_progress[-1, ]
  expect_warning(audpc_by_unit(short, "Severity_pct", "Day", c("Block", "Treatment")), "differ")
})

test_that("Smith's index: variances, fit and optimum plot size", {
  sm <- smith_index(uniformity_trial, cost_ratio = 2)
  expect_s3_class(sm, "ad_smith")
  A <- uniformity_trial
  tb <- sm$table
  expect_equal(tb$variance[tb$size == 1], stats::var(as.vector(A)))
  # 2 rows x 4 columns: plot means with explicit loops
  m <- c()
  for (i in seq(1, 24, 2)) for (j in seq(1, 16, 4)) m <- c(m, mean(A[i:(i + 1), j:(j + 3)]))
  expect_equal(tb$variance[tb$length == 2 & tb$width == 4], stats::var(m))
  u <- tb[tb$used, ]
  expect_true(all(u$length < 24 & u$width < 16))
  co <- stats::coef(stats::lm(log(variance) ~ log(size), u, weights = df))
  expect_equal(sm$b, unname(-co[2]))
  expect_equal(sm$V1, unname(exp(co[1])))
  expect_equal(sm$optimum, sm$b * 2 / (1 - sm$b))
  expect_equal(sum(smith_index(A, strips = TRUE)$table$used), nrow(tb))
  expect_output(AgriDesign:::print.ad_smith(sm), "Optimum plot size")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(nrow(AgriDesign:::plot.ad_smith(sm)), nrow(tb))
  expect_error(smith_index(matrix(1:4, 2)), "too few plot sizes")
  expect_error(smith_index(A[, 1:3] * NA), "missing")
})

test_that("Hatheway plot size follows its formula", {
  sm <- smith_index(uniformity_trial)
  h <- hatheway_size(sm, reps = c(3, 5), difference = 10, treatments = 6)
  t1 <- stats::qt(0.975, (h$reps - 1) * 5); t2 <- stats::qt(0.8, (h$reps - 1) * 5)
  expect_equal(h$plot_size, (2 * (t1 + t2)^2 * sm$cv1^2 / (h$reps * 100))^(1 / sm$b))
  expect_true(h$plot_size[2] < h$plot_size[1])
  hn <- hatheway_size(0.5, reps = 4, difference = 10, cv = 10)
  expect_equal(hn$plot_size, (2 * (stats::qnorm(0.975) + stats::qnorm(0.8))^2 * 100 / 400)^2)
  expect_error(hatheway_size(0.5, reps = 4, difference = 10), "coefficient of variation")
})
