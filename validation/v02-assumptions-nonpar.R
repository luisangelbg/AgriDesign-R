# Stage 2 validation: assumptions, transformations, Box-Cox and rank-based tests
# against base R, nortest, car, MASS, PMCMRplus, ARTool, rcompanion and the
# figures of the web application (user manual, chapter 5).
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages({ library(nortest); library(car); library(PMCMRplus); library(ARTool); library(rcompanion) })
ok <- 0; bad <- 0
chk <- function(what, got, ref, tol = 1e-6) {
  got <- as.numeric(got); ref <- as.numeric(ref)
  d <- max(abs(got - ref) / pmax(1, abs(ref)))
  pass <- length(got) == length(ref) && is.finite(d) && d < tol
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-58s got %-26s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(signif(got, 6), collapse = " "), paste(signif(ref, 6), collapse = " ")))
}
same <- function(what, got, ref) {
  pass <- identical(as.character(got), as.character(ref))
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-58s got %-26s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(got, collapse = " "), paste(ref, collapse = " ")))
}
tst <- function(a, name, col = "p") a$tests[a$tests$test == name, col]

# ---------------- maize (manual: grade A, SW W 0.930 p 0.2477, Levene median p 0.6672,
# Levene mean p 0.0287, Tukey p 0.8208, file row 17 stud 2.646, Box-Cox 0.70 [0, 1.5],
# Anderson-Darling A2 = 0.495, p = 0.1840)
f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
a <- check_assumptions(f)
e <- a$residuals$residual; g <- factor(maize_nitrogen$Nitrogen_kg_ha)
same("maize: grade A", a$grade, "A")
chk("maize: Shapiro W, p vs app", round(c(tst(a, "Shapiro-Wilk", "statistic"), tst(a, "Shapiro-Wilk")), 4), c(0.9300, 0.2477), 1e-3)
chk("maize: Levene median p vs app", round(tst(a, "Levene (median)"), 4), 0.6672)
chk("maize: Levene mean p vs app", round(tst(a, "Levene (mean)"), 4), 0.0287)
chk("maize: Levene median vs car::leveneTest", tst(a, "Levene (median)", "statistic"), car::leveneTest(e, g)[1, "F value"])
chk("maize: Levene mean vs car::leveneTest", tst(a, "Levene (mean)", "statistic"), car::leveneTest(e, g, center = mean)[1, "F value"])
chk("maize: Tukey additivity p vs app", round(tst(a, "Tukey one-degree-of-freedom"), 4), 0.8208)
chk("maize: Anderson-Darling vs nortest", c(tst(a, "Anderson-Darling", "statistic"), tst(a, "Anderson-Darling")),
    c(ad.test(e)$statistic, ad.test(e)$p.value), 1e-10)
chk("maize: A2 and p vs app", round(c(tst(a, "Anderson-Darling", "statistic"), tst(a, "Anderson-Darling")), 4), c(0.495, 0.1840), 1e-3)
lmm <- lm(Yield_t_ha ~ factor(Block) + g, maize_nitrogen)
chk("maize: studentized residuals vs rstandard", a$residuals$studentized, rstandard(lmm))
chk("maize: Durbin-Watson vs car", tst(a, "Durbin-Watson (data order)", "statistic"), car::durbinWatsonTest(lmm)$dw)
s <- a$suspicious
chk("maize: suspicious obs = data row 16 (file row 17), r = 2.646", c(s$row, round(s$studentized, 3)), c(16, 2.646))
bc <- box_cox(f)
chk("maize: Box-Cox lambda, interval vs app", c(bc$lambda, bc$lower, bc$upper), c(0.7, 0, 1.5))
mb <- MASS::boxcox(lmm, lambda = seq(-2, 2, 0.1), plotit = FALSE)
chk("maize: Box-Cox profile vs MASS (up to a constant)", diff(bc$profile$loglik), diff(mb$y), 1e-8)
# Tukey additivity: classic RCBD formula
ym <- tapply(maize_nitrogen$Yield_t_ha, list(maize_nitrogen$Block, g), mean)
bi <- rowMeans(ym) - mean(ym); tj <- colMeans(ym) - mean(ym)
ssn <- sum(ym * outer(bi, tj))^2 / (sum(bi^2) * sum(tj^2))
chk("maize: Tukey 1-df SS vs textbook formula", f$X |> (\(X) AgriDesign:::.tukey_additivity(X, f$y)[["ss"]])(), ssn)

# ---------------- tomato (manual: grade C by Tukey F 54.25, Box-Cox 0.20 [0.1, 0.4],
# suggested ln(y), grade A after using it)
f <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block")
a <- check_assumptions(f)
same("tomato: grade C", a$grade, "C")
chk("tomato: Tukey additivity F vs app", round(tst(a, "Tukey one-degree-of-freedom", "statistic"), 2), 54.25)
bc <- box_cox(f)
chk("tomato: Box-Cox vs app", c(bc$lambda, bc$lower, bc$upper), c(0.2, 0.1, 0.4))
tr <- transformations(f)
same("tomato: suggested transformation = log", tr$suggested, "log")
tl <- tomato_factorial; tl$lnY <- log(tl$Yield_kg_plant)
same("tomato: grade A after ln(y)", check_assumptions(design_anova(tl, "lnY", c("Variety", "Fertilizer"), block = "Block"))$grade, "A")
d <- tomato_factorial; grp <- interaction(d$Variety, d$Fertilizer)
chk("tomato: Bartlett vs bartlett.test", tst(a, "Bartlett", "statistic"), bartlett.test(d$Yield_kg_plant, grp)$statistic)
chk("tomato: Fligner vs fligner.test", tst(a, "Fligner-Killeen", "statistic"), fligner.test(d$Yield_kg_plant, grp)$statistic)

# ---------------- fungicides (manual: Friedman chi2 14.04, p 0.0029, W 0.936, letters a/b/c/d)
f <- design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block")
np <- nonparametric(f)
same("fungicides: method picked = friedman", np$method, "friedman")
ft <- friedman.test(Severity_score ~ Fungicide | Block, disease_scores)
chk("fungicides: Friedman chi2 vs friedman.test", np$test$statistic, ft$statistic, 1e-10)
chk("fungicides: chi2, p, W vs app", round(c(np$test$statistic, np$test$p, np$test$effect), c(2, 4, 3)), c(14.04, 0.0029, 0.936))
same("fungicides: four groups a/b/c/d", sort(np$groups$group), c("a", "b", "c", "d"))
ds <- disease_scores
pc <- frdAllPairsConoverTest(ds$Severity_score, factor(ds$Fungicide, levels = unique(ds$Fungicide)), factor(ds$Block), p.adjust.method = "none")
pm <- pc$p.value; lv <- levels(f$data$Fungicide)
ref <- mapply(function(i, j) { v <- pm[i, j]; if (is.na(v)) pm[j, i] else v }, np$pairs$level_j, np$pairs$level_i)
chk("fungicides: Conover p vs PMCMRplus", np$pairs$p, ref, 1e-8)
nm <- nonparametric(f, posthoc = "nemenyi")
pn <- frdAllPairsNemenyiTest(ds$Severity_score, factor(ds$Fungicide, levels = unique(ds$Fungicide)), factor(ds$Block))$p.value
ref <- mapply(function(i, j) { v <- pn[i, j]; if (is.na(v)) pn[j, i] else v }, nm$pairs$level_j, nm$pairs$level_i)
chk("fungicides: Nemenyi p vs PMCMRplus", nm$pairs$p, ref, 1e-6)

# ---------------- beans (manual: Kruskal-Wallis H 17.63, p 0.0015)
f <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
np <- nonparametric(f)
kt <- kruskal.test(Yield_kg_plot ~ Variety, bean_varieties)
chk("beans: KW H vs kruskal.test", np$test$statistic, kt$statistic, 1e-10)
chk("beans: KW H, p vs app", round(c(np$test$statistic, np$test$p), c(2, 4)), c(17.63, 0.0015))
bv <- bean_varieties; bv$Variety <- factor(bv$Variety, levels = levels(f$data$Variety))
pd <- kwAllPairsDunnTest(bv$Yield_kg_plot, bv$Variety, p.adjust.method = "none")$p.value
ref <- mapply(function(i, j) { v <- pd[i, j]; if (is.na(v)) pd[j, i] else v }, np$pairs$level_j, np$pairs$level_i)
chk("beans: Dunn p vs PMCMRplus", np$pairs$p, ref, 1e-8)
chk("beans: Holm adjustment vs p.adjust", np$pairs$p_adj, p.adjust(np$pairs$p, "holm"))
w <- nonparametric(f, method = "welch")
ow <- oneway.test(Yield_kg_plot ~ Variety, bean_varieties)
chk("beans: Welch F, df2 vs oneway.test", c(w$test$statistic, w$test$df[2]), c(ow$statistic, ow$parameter[2]))
gh <- gamesHowellTest(bv$Yield_kg_plot, bv$Variety)$p.value
ref <- mapply(function(i, j) { v <- gh[i, j]; if (is.na(v)) gh[j, i] else v }, w$pairs$level_j, w$pairs$level_i)
chk("beans: Games-Howell p vs PMCMRplus", w$pairs$p, ref, 1e-6)
mw <- nonparametric(f, method = "kruskal_mw")
w1 <- suppressWarnings(wilcox.test(bv$Yield_kg_plot[bv$Variety == mw$pairs$level_i[1]],
                                   bv$Yield_kg_plot[bv$Variety == mw$pairs$level_j[1]], exact = FALSE, correct = FALSE))
chk("beans: Mann-Whitney p (first pair) vs wilcox.test", mw$pairs$p[1], w1$p.value)

# ---------------- factorial in CRD: ART vs ARTool, SRH vs rcompanion
set.seed(11)
fd <- expand.grid(rep = 1:4, A = c("a1", "a2", "a3"), B = c("b1", "b2"))
fd$y <- round(rexp(nrow(fd), 1 / (3 + as.integer(fd$A) + 2 * (fd$B == "b2"))), 2)
f <- design_anova(fd, "y", c("A", "B"))
art_np <- nonparametric(f, method = "art")
m <- art(y ~ A * B, data = fd)
aa <- anova(m)
chk("ART: F of A, B, A:B vs ARTool", art_np$table$F, aa$F, 1e-8)
srh <- nonparametric(f, method = "srh")
sr <- scheirerRayHare(y ~ A + B, data = fd, verbose = FALSE)
chk("SRH: H of A, B, A:B vs rcompanion", srh$table$H, sr$H[1:3], 1e-8)

# ---------------- compact letters: a non-transitive case
dm <- matrix(FALSE, 3, 3); dm[1, 3] <- dm[3, 1] <- TRUE
same("letters: a, ab, b when only 1 and 3 differ", AgriDesign:::.letters(c(3, 2, 1), dm), c("a", "ab", "b"))

cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
