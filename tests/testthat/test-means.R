fac <- function(d, v) { for (x in v) d[[x]] <- factor(d[[x]]); d }

test_that("Tukey reproduces TukeyHSD and the figures of the web application", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  tk <- compare_means(f, "Nitrogen_kg_ha")
  d <- fac(maize_nitrogen, c("Block", "Nitrogen_kg_ha"))
  th <- stats::TukeyHSD(stats::aov(Yield_t_ha ~ Block + Nitrogen_kg_ha, d), "Nitrogen_kg_ha")[[1]]
  expect_equal(tk$pairs$p, unname(th[paste(tk$pairs$level_j, tk$pairs$level_i, sep = "-"), "p adj"]),
               tolerance = 1e-6)
  expect_equal(round(unname(tk$msd), 3), 0.283)
  expect_equal(tk$means$group, c("d", "c", "b", "a"))
})

test_that("multiple range tests give the critical ranges of the textbook formulas", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  se <- sqrt(f$residual[["ms"]] / 4); df <- f$residual[["df"]]
  d <- compare_means(f, "Nitrogen_kg_ha", method = "duncan")
  expect_equal(unname(d$msd), stats::qtukey(0.95^(1:3), 2:4, df) * se)
  s <- compare_means(f, "Nitrogen_kg_ha", method = "snk")
  expect_equal(unname(s$msd), stats::qtukey(0.95, 2:4, df) * se)
  r <- compare_means(f, "Nitrogen_kg_ha", method = "regwq")
  expect_equal(unname(r$msd), stats::qtukey(1 - c(1 - 0.95^(2 / 4), 0.05, 0.05), 2:4, df) * se)
  sc <- compare_means(f, "Nitrogen_kg_ha", method = "scheffe")
  expect_equal(unname(sc$msd), sqrt(3 * stats::qf(0.95, 3, df)) * sqrt(2) * se)
  l <- compare_means(f, "Nitrogen_kg_ha", method = "lsd")
  expect_equal(unname(l$msd), stats::qt(0.975, df) * sqrt(2) * se)
})

test_that("Bonferroni, Holm and Sidak follow pairwise.t.test", {
  f <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
  bv <- bean_varieties; bv$V <- factor(bv$Variety, levels = levels(f$data$Variety))
  pm <- stats::pairwise.t.test(bv$Yield_kg_plot, bv$V, p.adjust.method = "none")$p.value
  get <- function(a, b) mapply(function(x, y) { v <- pm[x, y]; if (is.na(v)) pm[y, x] else v }, a, b,
                               USE.NAMES = FALSE)
  for (m in c("bonferroni", "holm")) {
    cm <- compare_means(f, "Variety", method = m)
    expect_equal(cm$pairs$p_adj, stats::p.adjust(get(cm$pairs$level_j, cm$pairs$level_i), m))
  }
  si <- compare_means(f, "Variety", method = "sidak")
  expect_equal(si$pairs$p_adj, 1 - (1 - get(si$pairs$level_j, si$pairs$level_i))^10)
  expect_equal(round(unname(compare_means(f, "Variety", method = "bonferroni")$msd), 3), 0.239)
})

test_that("Dunnett's critical value matches the web application and the tables", {
  f <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
  dn <- compare_means(f, "Variety", method = "dunnett", control = "Negro Jamapa")
  expect_equal(round(unname(dn$msd), 3), 0.198)
  expect_equal(round(dn$pairs$p[dn$pairs$level_i == "Bayo"], 4), 0.0318)
  expect_equal(dn$means$group[dn$means$label == "Negro Jamapa"], "control")
  # equicorrelated rho = 0.5, 4 comparisons, df = 15: two-sided 5 % point 2.7273
  expect_equal(AgriDesign:::.dunnett_crit(0.05, rep(sqrt(0.5), 4), 15), 2.7273, tolerance = 1e-4)
  # one comparison reduces to Student's t
  expect_equal(AgriDesign:::.dunnett_p(2.1, sqrt(0.5), 10), 2 * stats::pt(-2.1, 10), tolerance = 1e-8)
  expect_error(compare_means(f, "Variety", method = "dunnett", control = "Pinto"), "not found")
})

test_that("split-plot comparisons use the right errors and Satterthwaite df", {
  f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "split_rcbd")
  ea <- f$strata[["Block:Irrigation"]]; eb <- f$strata$Residuals
  mi <- compare_means(f, "Irrigation", method = "lsd")
  expect_equal(mi$error$ms, ea$ms)
  expect_equal(unname(mi$msd), stats::qt(0.975, ea$df) * sqrt(2 * ea$ms / 12))
  w <- compare_means(f, "Variety", by = "Irrigation")
  expect_equal(w$error$ms, eb$ms)
  a <- compare_means(f, "Irrigation", by = "Variety")
  ms <- (3 * eb$ms + ea$ms) / 4
  df <- ms^2 / ((3 * eb$ms / 4)^2 / eb$df + (ea$ms / 4)^2 / ea$df)
  expect_equal(a$error$ms, ms)
  expect_equal(round(a$error$df, 1), 10.7)
  expect_error(compare_means(f, c("Irrigation", "Variety")), "simple effects")
})

test_that("trend and contrasts reproduce the web application", {
  f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
  tr <- trend(f, "Nitrogen_kg_ha")
  expect_equal(round(tr$components$ss[1:2], 4), c(19.7707, 1.7756))
  expect_equal(tr$components$coefficients[1], "-3, -1, 1, 3")
  expect_equal(round(tr$optimum[["x"]], 1), 179.5)
  expect_equal(sum(tr$components$ss), f$table$ss[f$table$term == "Nitrogen_kg_ha"])
  ct <- test_contrasts(f, "Nitrogen_kg_ha", c(0, 0, 1, -1))
  expect_equal(round(ct$p, 4), 0.0095)
  expect_error(test_contrasts(f, "Nitrogen_kg_ha", c(1, 1, 0, 0)), "sum to zero")
  # unequally spaced doses: polynomial coefficients stay orthogonal
  P <- AgriDesign:::.ortho_poly(c(0, 50, 100, 200), rep(4, 4), 3)
  expect_equal(sum(P[[1]] * P[[2]]), 0)
  expect_equal(sum(P[[1]]), 0)
})

test_that("transformed analyses report back-transformed means", {
  f <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
                    block = "Block", transform = "log")
  m <- ls_means(f, "Variety")
  expect_equal(m$back_transformed, exp(m$mean))
  expect_error(design_anova(tomato_factorial, "Yield_kg_plant", "Variety", transform = "logit"),
               "does not suit")
})
