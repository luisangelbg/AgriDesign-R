# Plant-breeding analyses: AMMI, stability, line x tester, diallels, North Carolina designs.
# References are base R fits and the closed formulas of the original papers.
# Print and plot methods are called directly, so the tests do not depend on S3 registration.

row_of <- function(tab, s) tab[trimws(tab$source) == s, ]

test_that("AMMI reproduces lm() and the singular value decomposition of the interaction", {
  fit <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
  d <- maize_met; for (v in c("Hybrid", "Environment", "Block")) d[[v]] <- factor(d[[v]])
  a <- stats::anova(stats::lm(stats::terms(Yield_t_ha ~ Environment + Environment:Block + Hybrid + Hybrid:Environment,
                                           keep.order = TRUE), d))
  expect_equal(fit$anova$ss[1:4], a[["Sum Sq"]][1:4])
  expect_equal(fit$mse, a["Residuals", "Mean Sq"])
  expect_equal(row_of(fit$anova, "ENV")$F, a[1, "Mean Sq"] / a[2, "Mean Sq"])
  Y <- fit$means
  Z <- Y - outer(rowMeans(Y), colMeans(Y), "+") + mean(Y)
  expect_equal(fit$axes$singular_value, svd(Z)$d[1:5])
  G <- as.matrix(fit$genotypes[paste0("IPCA", 1:5)]); E <- as.matrix(fit$environments[paste0("IPCA", 1:5)])
  expect_equal(unname(G %*% t(E)), unname(Z))
  expect_equal(sum(fit$axes$ss), row_of(fit$anova, "GxE")$ss)
  expect_equal(sum(fit$axes$df), 11 * 5)
  ss <- fit$axes$ss
  expect_equal(fit$genotypes$asv, sqrt((ss[1] / ss[2] * G[, 1])^2 + G[, 2]^2))
  expect_equal(fit$genotypes$ysi, rank(fit$genotypes$asv) + rank(-fit$genotypes$mean))
  f2 <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block", n_axes = 2)
  expect_equal(row_of(f2$anova, "Residual GxE")$ss, sum(ss[3:5]))
})

test_that("AMMI on means with the error of the trial equals the plot analysis", {
  fb <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
  mm <- stats::aggregate(Yield_t_ha ~ Hybrid + Environment, maize_met, mean)
  fm <- ammi(mm, "Yield_t_ha", "Hybrid", "Environment", mse = fb$mse, df_error = fb$df_error, reps = 3)
  expect_equal(fm$axes[c("ss", "F", "p")], fb$axes[c("ss", "F", "p")])
  expect_equal(row_of(fm$anova, "GEN")$ss, row_of(fb$anova, "GEN")$ss)
  expect_true(all(is.na(ammi(mm, "Yield_t_ha", "Hybrid", "Environment")$axes$F)))
})

test_that("AMMI checks its input", {
  expect_error(ammi(maize_met[-1, ], "Yield_t_ha", "Hybrid", "Environment"), "same in every")
  d <- maize_met[!(maize_met$Hybrid == "H01" & maize_met$Environment == "E1"), ]
  expect_error(ammi(d, "Yield_t_ha", "Hybrid", "Environment"), "missing: H01 x E1")
  expect_error(ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", mse = 1), "one value")
  expect_error(ammi(maize_met, "Yield", "Hybrid", "Environment"), "not found")
  expect_error(ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", n_axes = 9), "n_axes")
})

test_that("AMMI biplots draw and restore the graphical parameters", {
  fit <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  op <- graphics::par("mar")
  p1 <- plot.ad_ammi(fit)
  p2 <- plot.ad_ammi(fit, type = "ammi2")
  expect_equal(graphics::par("mar"), op)
  expect_equal(nrow(p1), 12 + 6)
  expect_equal(p2$x[p2$type == "environment"], fit$environments$IPCA1)
  expect_error(plot.ad_ammi(fit, type = "gge"))
})

test_that("Eberhart-Russell and Finlay-Wilkinson follow lm() of every genotype", {
  st <- stability(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
  Y <- st$means; I <- colMeans(Y) - mean(Y)
  fits <- lapply(seq_len(nrow(Y)), function(i) summary(stats::lm(Y[i, ] ~ I)))
  expect_equal(st$eberhart$table$b, vapply(fits, function(f) f$coefficients[2, 1], 0))
  expect_equal(st$eberhart$table$se_b, vapply(fits, function(f) f$coefficients[2, 2], 0))
  expect_equal(st$eberhart$table$s2d, vapply(fits, function(f) f$sigma^2, 0) - st$mse / 3)
  expect_equal(mean(st$eberhart$table$b), 1)
  x <- colMeans(Y)
  expect_equal(st$finlay$intercept, vapply(seq_len(nrow(Y)), function(i) unname(stats::coef(stats::lm(Y[i, ] ~ x))[1]), 0))
  tab <- st$anova
  expect_equal(row_of(tab, "Heterogeneity")$ss + row_of(tab, "Residual")$ss, row_of(tab, "GxE")$ss)
  expect_equal(row_of(tab, "Residual")$df, 11 * 4)
})

test_that("Shukla, Wricke and Kang statistics follow their definitions", {
  st <- stability(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
  Y <- st$means; g <- nrow(Y); e <- ncol(Y)
  Z <- Y - outer(rowMeans(Y), colMeans(Y), "+") + mean(Y)
  W <- 3 * rowSums(Z^2)
  expect_equal(st$wricke$W, unname(W))
  expect_equal(sum(st$wricke$W), row_of(st$anova, "GxE")$ss)
  expect_equal(st$shukla$sigma2, unname(g / ((g - 2) * (e - 1)) * W - sum(W) / ((g - 1) * (g - 2) * (e - 1))))
  expect_equal(sum(st$shukla$sigma2), g * row_of(st$anova, "GxE")$ms)   # sum of sigma2_i = g MS(GxE)
  k <- st$kang$table
  expect_equal(st$kang$lsd, stats::qt(0.975, st$df_error) * sqrt(2 * st$mse / (3 * e)))
  expect_equal(k$ys, k$rank + k$adjustment + k$rating)
  expect_equal(k$selected, k$ys > mean(k$ys))
  expect_equal(k$rating, ifelse(st$shukla$p < 0.01, -8, ifelse(st$shukla$p < 0.05, -4, ifelse(st$shukla$p < 0.1, -2, 0))))
  mm <- stats::aggregate(Yield_t_ha ~ Hybrid + Environment, maize_met, mean)
  expect_null(stability(mm, "Yield_t_ha", "Hybrid", "Environment", method = "kang")$kang)
})

test_that("Huehn's expectations and variances match the exact null distribution", {
  K <- 4; N <- 3
  R <- as.matrix(expand.grid(rep(list(seq_len(K)), N)))
  s1 <- apply(R, 1, function(v) sum(abs(outer(v, v, "-"))) / (N * (N - 1)))
  s2 <- apply(R, 1, stats::var)
  d <- data.frame(G = rep(seq_len(K), N), E = rep(seq_len(N), each = K), y = c(outer(seq_len(K), seq_len(N))) + 0.5)
  h <- stability(d, "y", "G", "E", method = "huehn")$huehn
  expect_equal(unname(h$expectation), c(mean(s1), mean((s1 - mean(s1))^2), mean(s2), mean((s2 - mean(s2))^2)))
  st <- stability(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block", method = "huehn")$huehn
  Y <- tapply(maize_met$Yield_t_ha, list(maize_met$Hybrid, maize_met$Environment), mean)
  Rk <- apply(Y - rowMeans(Y), 2, rank)
  expect_equal(st$table$S2, unname(apply(Rk, 1, stats::var)))
  expect_equal(unname(st$sum_z[1]), sum(st$table$Z1))
})

test_that("line x tester partitions the crosses as lm() does", {
  lt <- line_tester(maize_lxt, "Yield_t_ha", "Line", "Tester", block = "Block")
  d <- maize_lxt; d$Block <- factor(d$Block)
  d$entry <- factor(paste(d$Line, d$Tester))
  m <- stats::lm(Yield_t_ha ~ Block + entry, d)
  expect_equal(lt$mse, stats::deviance(m) / stats::df.residual(m))
  cr <- d[!is.na(d$Line) & !is.na(d$Tester), ]
  a <- stats::anova(stats::lm(Yield_t_ha ~ Block + Line * Tester, cr))
  expect_equal(c(row_of(lt$anova, "Lines")$ss, row_of(lt$anova, "Testers")$ss, row_of(lt$anova, "Lines x Testers")$ss),
               a[2:4, "Sum Sq"])
  tab <- lt$anova
  expect_equal(row_of(tab, "Parents")$ss + row_of(tab, "Parents vs crosses")$ss + row_of(tab, "Crosses")$ss,
               row_of(tab, "Entries")$ss)
  expect_equal(sum(lt$gca_lines$gca), 0); expect_equal(sum(lt$gca_testers$gca), 0)
  expect_equal(as.vector(tapply(lt$sca$sca, lt$sca$line, sum)), rep(0, 6))
  expect_equal(sum(lt$contribution), 100)
  expect_equal(lt$se[["gca_line"]], sqrt(lt$mse / (3 * 3)))
  cv <- lt$components$estimate
  expect_equal(cv[3], (5 * 3 * cv[1] + 2 * 6 * cv[2]) / (5 * 3 + 2 * 6))
  expect_equal(lt$genetic$sigma2_D, c(4, 1) * cv[6])
  expect_error(line_tester(maize_lxt[-1, ], "Yield_t_ha", "Line", "Tester", block = "Block"), "every block")
  expect_error(line_tester(subset(maize_lxt, Tester != "T3" | is.na(Tester)), "Yield_t_ha", "Line", "Tester",
                           block = "Block")$se, NA)
})

test_that("diallel methods 2 and 4 reproduce Griffing's closed formulas", {
  dd <- maize_diallel; P <- paste0("P", 1:6); p <- 6
  x <- tapply(dd$Yield_t_ha, list(factor(dd$Female, P), factor(dd$Male, P)), mean)
  d4 <- diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", method = 4)
  xu <- x; xu[!upper.tri(xu)] <- 0; xs <- xu + t(xu); Xi <- rowSums(xs); X <- sum(xu)
  expect_equal(d4$anova$ss[1:2], c(sum(Xi^2) / (p - 2) - 4 * X^2 / (p * (p - 2)),
                                   sum(xu^2) - sum(Xi^2) / (p - 2) + 2 * X^2 / ((p - 1) * (p - 2))))
  expect_equal(d4$gca$estimate, unname((p * Xi - 2 * X) / (p * (p - 2))))
  s2 <- d4$mse / 3
  expect_equal(d4$gca$se, rep(sqrt((p - 1) / (p * (p - 2)) * s2), p))
  expect_equal(d4$sca$se[1], sqrt((p - 3) / (p - 1) * s2))
  expect_equal(d4$se_differences$se, sqrt(c(2 / (p - 2), 2 * (p - 3) / (p - 2), 2 * (p - 4) / (p - 2)) * s2))
  expect_equal(d4$rows_left_out, (6 + 15) * 3)   # parents and reciprocals, three blocks
  d2 <- diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", method = 2)
  xu <- x; xu[lower.tri(xu)] <- 0; xs <- xu + t(xu); diag(xs) <- diag(x)
  Xi <- rowSums(xs); X <- sum(xu); xii <- diag(x)
  expect_equal(d2$gca$estimate, unname((Xi + xii - 2 * X / p) / (p + 2)))
  expect_equal(d2$anova$ss[1], (sum((Xi + xii)^2) - 4 * X^2 / p) / (p + 2))
  ms <- stats::setNames(d2$anova$ms, trimws(d2$anova$source))
  expect_equal(d2$model2$estimate[1:2], c((ms[["GCA"]] - ms[["SCA"]]) / (p + 2), ms[["SCA"]] - ms[["Error"]]))
})

test_that("diallel methods 1 and 3 split the reciprocal effects", {
  dd <- maize_diallel; P <- paste0("P", 1:6); p <- 6
  x <- tapply(dd$Yield_t_ha, list(factor(dd$Female, P), factor(dd$Male, P)), mean)
  d1 <- diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", method = 1)
  ss <- stats::setNames(d1$anova$ss, trimws(d1$anova$source))
  expect_equal(ss[["Reciprocal"]], sum(((x - t(x))[upper.tri(x)])^2) / 2)
  expect_equal(ss[["Maternal"]] + ss[["Non-maternal"]], ss[["Reciprocal"]])
  Xi <- rowSums(x) + colSums(x)
  expect_equal(d1$gca$estimate, unname(Xi / (2 * p) - sum(x) / p^2))
  expect_equal(d1$reciprocal$estimate[1], (x[1, 2] - x[2, 1]) / 2)
  xo <- x; diag(xo) <- 0
  expect_equal(d1$maternal$estimate, unname((rowSums(xo) - colSums(xo)) / (2 * p)))
  d <- dd; d$B <- factor(d$Block); d$E <- factor(paste(d$Female, d$Male))
  a <- stats::anova(stats::lm(Yield_t_ha ~ B + E, d))
  expect_equal(d1$mse, a["Residuals", "Mean Sq"])
  d3 <- diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", method = 3)
  expect_equal(d3$anova$df[1:3], c(5, 9, 15))
  expect_equal(d3$rows_left_out, 18)
})

test_that("diallel checks its input", {
  dd <- maize_diallel
  expect_error(diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", method = 5), "1, 2, 3 or 4")
  bad <- dd[!(dd$Female == "P1" & dd$Male == "P2"), ]
  expect_error(diallel(bad, "Yield_t_ha", "Female", "Male", block = "Block", method = 3), "P1 x P2")
  small <- dd[dd$Female %in% c("P1", "P2", "P3") & dd$Male %in% c("P1", "P2", "P3"), ]
  expect_error(diallel(small, "Yield_t_ha", "Female", "Male", block = "Block", method = 4), "four parents")
  expect_error(diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", inbreeding = 2), "inbreeding")
})

test_that("North Carolina designs follow lm() and their expected mean squares", {
  d1 <- nc_designs$nc1
  f1 <- nc_design(d1, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 1)
  d <- d1; d$S <- factor(d$Set); d$R <- factor(paste(d$Set, d$Rep)); d$M <- factor(paste(d$Set, d$Male))
  d$F <- factor(paste(d$Set, d$Male, d$Female))
  a <- stats::anova(stats::lm(Yield_t_ha ~ S + R + M + F, d))
  expect_equal(f1$anova$ss[1:5], a[["Sum Sq"]])
  ms <- a[["Mean Sq"]]
  expect_equal(f1$components$estimate, c((ms[3] - ms[4]) / (2 * 4), (ms[4] - ms[5]) / 2))
  expect_equal(row_of(f1$anova, "Males(sets)")$F, ms[3] / ms[4])
  f2 <- nc_design(nc_designs$nc2, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 2)
  ms2 <- stats::setNames(f2$anova$ms, f2$anova$source)
  expect_equal(f2$components$estimate[3], (ms2[["Males x Females(sets)"]] - ms2[["Error"]]) / 2)
  expect_equal(f2$genetic$estimate[3], 4 * mean(f2$components$estimate[1:2]))
  d3 <- nc_designs$nc3
  f3 <- nc_design(d3, "Yield_t_ha", "Male", "Tester", block = "Rep", set = "Set", design = 3)
  pm <- tapply(d3$Yield_t_ha, list(paste(d3$Set, d3$Male), d3$Tester), mean)
  setm <- substr(rownames(pm), 1, 1); sm <- pm[, 1] + pm[, 2]; df_ <- pm[, 1] - pm[, 2]
  vS <- sum(tapply(sm, setm, function(v) sum((v - mean(v))^2))) / (12 - 2)
  vD <- sum(tapply(df_, setm, function(v) sum((v - mean(v))^2))) / (12 - 2)
  mse <- row_of(f3$anova, "Error")$ms
  expect_equal(f3$genetic$estimate[1:2], c(vS - mse, (vD - mse) / 2))
  expect_error(nc_design(d3, "Yield_t_ha", "Male", "Tester", block = "Rep", design = 4), "1, 2 or 3")
  d3b <- d3; d3b$Tester[1] <- "P3"
  expect_error(nc_design(d3b, "Yield_t_ha", "Male", "Tester", block = "Rep", set = "Set", design = 3), "two testers")
})

test_that("North Carolina design I separates the plot error from plants within plots", {
  d <- nc_designs$nc1
  d2 <- rbind(d, transform(d, Yield_t_ha = Yield_t_ha + rep(c(-0.1, 0.1), length.out = nrow(d))))
  f <- nc_design(d2, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 1)
  expect_equal(unname(f$sizes["per_plot"]), 2)
  expect_true(all(c("Plot error", "Within plots") %in% f$anova$source))
  ms <- stats::setNames(f$anova$ms, f$anova$source)
  expect_equal(f$components$estimate[2], (ms[["Females(males)"]] - ms[["Plot error"]]) / (2 * 2))
})

test_that("print methods run", {
  expect_output(print.ad_ammi(ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")), "IPCA1")
  expect_output(print.ad_stability(stability(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")), "Kang")
  expect_output(print.ad_line_tester(line_tester(maize_lxt, "Yield_t_ha", "Line", "Tester", block = "Block")), "GCA")
  expect_output(print.ad_diallel(diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block")), "Model II")
  expect_output(print.ad_nc_design(nc_design(nc_designs$nc2, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set",
                                           design = 2)), "Expected mean squares")
})
