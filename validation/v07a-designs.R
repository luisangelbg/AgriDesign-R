# Stage 7a validation: Graeco-Latin, Youden, BIB, cyclic, alpha and lattice designs
# (structure over many seeds, efficiency factors against closed forms and against
# agricolae's designs), their analysis against aov()/agricolae, and the REML
# recovery of inter-block information against lme4 + emmeans.
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages({ library(agricolae); library(lme4); library(emmeans) })
ok <- 0; bad <- 0
chk <- function(what, pass, info = "") {
  pass <- isTRUE(all(pass)); if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %s %s\n", if (pass) "OK" else "FAIL", what, info))
}
near <- function(a, b, tol = 1e-6) length(a) == length(b) && max(abs(as.numeric(a) - as.numeric(b)) / pmax(1, abs(as.numeric(b)))) < tol
all_seeds <- function(n, f) all(vapply(seq_len(n), function(s) isTRUE(all(f(s))), NA))
concurrence <- function(trt, blk) { N <- unclass(table(trt, blk)); M <- N %*% t(N); M[upper.tri(M)] }

# ---------------- structure ----------------
for (t in c(3, 4, 5, 7, 8, 9)) chk(sprintf("Graeco-Latin %d x %d: Latin in both letters and orthogonal (100 seeds)", t, t), all_seeds(100, function(s) {
  b <- make_design("graeco", trt = t, seed = s)$book
  all(table(b$Row, b$Treatment) == 1) && all(table(b$Column, b$Treatment) == 1) &&
    all(table(b$Row, b$Greek) == 1) && all(table(b$Column, b$Greek) == 1) && all(table(b$Treatment, b$Greek) == 1) }))
chk("Graeco-Latin of order 6 and 2 refused", inherits(try(make_design("graeco", trt = 6, seed = 1), silent = TRUE), "try-error") &&
      inherits(try(make_design("graeco", trt = 2, seed = 1), silent = TRUE), "try-error"))
for (key in names(AgriDesign:::.difference_sets)) {
  v <- as.integer(sub("-.*", "", key)); D <- AgriDesign:::.difference_sets[[key]]; k <- length(D)
  dif <- outer(D, D, "-") %% v; dif <- dif[row(dif) != col(dif)]
  lam <- k * (k - 1) / (v - 1)
  chk(sprintf("Difference set (%d, %d, %g): every non-zero difference %g times", v, k, lam, lam),
      lam == round(lam) && all(tabulate(dif, v - 1) == lam))
}
for (tk in list(c(7, 3), c(7, 4), c(11, 5), c(13, 4), c(6, 5))) chk(sprintf("Youden t = %d, k = %d: rows complete, columns balanced (50 seeds)", tk[1], tk[2]), all_seeds(50, function(s) {
  b <- make_design("youden", trt = tk[1], k = tk[2], seed = s)$book
  cc <- concurrence(b$Treatment, b$Column)
  all(table(b$Row, b$Treatment) == 1) && all(table(b$Column) == tk[2]) && length(unique(cc)) == 1 }))
for (tk in list(c(7, 3), c(6, 3), c(5, 2), c(9, 4))) chk(sprintf("BIB t = %d, k = %d: equal replication and equal concurrence (50 seeds)", tk[1], tk[2]), all_seeds(50, function(s) {
  b <- make_design("bib", trt = tk[1], k = tk[2], seed = s)$book
  length(unique(table(b$Treatment))) == 1 && length(unique(concurrence(b$Treatment, b$Block))) == 1 &&
    all(table(b$Block) == tk[2]) && all(table(b$Block, b$Treatment) <= 1) }))
chk("BIB efficiency = lambda t / (r k)", {
  lay <- make_design("bib", trt = 7, k = 3, seed = 1); near(design_efficiency(lay)$E, 1 * 7 / (3 * 3)) })
chk("Complete blocks have efficiency 1", near(design_efficiency(make_design("rcbd", trt = 6, r = 3, seed = 1))$E, 1))
for (k in c(3, 5)) for (r in 2:3) {
  lay <- make_design("lattice", trt = k^2, k = k, r = r, seed = 1)
  E <- (k^2 - 1) / (r * (k - 1) * r / (r - 1) + (k - 1) * (k + 1 - r))
  chk(sprintf("Square lattice %d x %d, r = %d: concurrences 0/1 and E = closed form %.4f", k, k, r, E),
      all(concurrence(lay$book$Treatment, paste(lay$book$Block, lay$book$IncBlock)) <= 1) && near(design_efficiency(lay)$E, E))
}
chk("Alpha designs resolvable with blocks of k (50 seeds)", all_seeds(50, function(s) {
  b <- make_design("alpha", trt = 20, k = 4, r = 3, seed = s)$book
  all(table(b$Block, b$Treatment) == 1) && all(table(paste(b$Block, b$IncBlock)) == 4) }))

# ---------------- efficiency against agricolae's designs ----------------
cat("\nEfficiency factor: our computation of agricolae's own layouts, and our designs\n")
for (p in list(c(12, 3, 2), c(20, 4, 2), c(20, 4, 3), c(30, 5, 3), c(24, 4, 2), c(35, 5, 2))) {
  t <- p[1]; k <- p[2]; r <- p[3]
  ag <- suppressWarnings(utils::capture.output(ag <- design.alpha(paste0("G", 1:t), k = k, r = r, seed = 1)))
  ag <- design.alpha(paste0("G", 1:t), k = k, r = r, seed = 1)
  bk <- ag$book
  E_true <- .efficiency(bk[[4]], paste(bk$replication, bk$block))$E   # column 4: the treatments
  E_printed <- ag$statistics$Efficiency
  lay <- make_design("alpha", trt = t, k = k, r = r, seed = 1); e <- design_efficiency(lay)
  # agricolae reports the Patterson-Williams upper bound as "Efficiency factor"; the layout itself is evaluated here
  chk(sprintf("alpha t=%d k=%d r=%d: agricolae reports the upper bound (%.4f); its layout has E = %.4f", t, k, r, E_printed, E_true),
      near(E_printed, e$upper_bound, 1e-6) && E_true <= E_printed + 1e-9)
  chk(sprintf("alpha t=%d k=%d r=%d: ours %.4f >= agricolae's layout %.4f, <= bound %.4f", t, k, r, e$E, E_true, e$upper_bound),
      e$E >= E_true - 1e-9 && e$E <= e$upper_bound + 1e-9)
}
for (p in list(c(10, 3, 6), c(12, 4, 8), c(9, 3, 6))) {
  t <- p[1]; k <- p[2]; r <- p[3]
  ag <- suppressWarnings(design.cyclic(paste0("G", 1:t), k = k, r = r, seed = 1))
  bk <- ag$book                                       # plots, group (initial block), block, treatment
  E_ag <- .efficiency(bk[[4]], paste(bk$group, bk$block))$E
  E_us <- design_efficiency(make_design("cyclic", trt = t, k = k, r = r, seed = 1))$E
  chk(sprintf("cyclic t=%d k=%d r=%d: ours %.4f >= agricolae's layout %.4f", t, k, r, E_us, E_ag), E_us >= E_ag - 1e-9)
}

# ---------------- analysis ----------------
sim <- function(lay, sb = 1) {
  b <- field_book(lay, "y"); set.seed(5)
  blk <- if ("IncBlock" %in% names(b)) paste(b$Block, b$IncBlock) else if ("Block" %in% names(b)) b$Block else b$Column
  b$y <- 10 + as.integer(factor(b$Treatment)) / 4 + rnorm(length(unique(blk)), 0, sb)[as.integer(factor(blk))] + rnorm(nrow(b), 0, 0.6)
  b
}
b <- sim(make_design("graeco", trt = 5, seed = 3)); b$Greek <- factor(b$Greek)
f <- design_anova(b, "y", "Treatment", block = "Greek", row = "Row", col = "Column")
a <- summary(aov(y ~ factor(Row) + factor(Column) + Greek + Treatment, b))[[1]]
chk("Graeco-Latin: suggested design and F of treatments vs aov", f$design == "graeco" && near(f$table$F[4], a[4, "F value"]))
b <- sim(make_design("youden", trt = 7, k = 3, seed = 3))
f <- design_anova(b, "y", "Treatment", row = "Row", col = "Column")
l1 <- lm(y ~ factor(Row) + factor(Column), b); l2 <- lm(y ~ factor(Row) + factor(Column) + Treatment, b)
chk("Youden: suggested design and adjusted SS of treatments vs lm", f$design == "youden" &&
      near(f$table$ss[f$table$term == "Treatment"], anova(l1, l2)[2, "Sum of Sq"]))
b <- sim(make_design("bib", trt = 7, k = 3, seed = 3))
f <- design_anova(b, "y", "Treatment", block = "Block")
bt <- with(b, BIB.test(factor(Block), Treatment, y, test = "tukey", console = FALSE))
lm1 <- lm(y ~ factor(Block) + Treatment, b, contrasts = list(Treatment = contr.sum))
em <- summary(emmeans(lm1, "Treatment"))
chk("BIB: suggested design, adjusted means vs emmeans", f$design == "bib" &&
      near(ls_means(f, "Treatment")$mean, em$emmean[match(levels(f$data$Treatment), em$Treatment)]))
chk("BIB: efficiency and CV vs agricolae BIB.test",
      near(design_efficiency(b, "Treatment", "Block")$E, bt$statistics$Efficiency) && near(f$cv, bt$statistics$CV, 1e-4))
chk("BIB: adjusted means vs agricolae BIB.test",
      near(ls_means(f, "Treatment")$mean, bt$means[levels(f$data$Treatment), "mean.adj"]))
# REML recovery: alpha design, lattice, BIBD and augmented against lme4
rec_chk <- function(label, b, rand, ...) {
  f <- design_anova(b, "y", "Treatment", ...)
  rc <- recover_interblock(f)
  m <- lmer(as.formula(paste("y ~", if (length(list(...)$block) == 2) "factor(Block) + " else "", "Treatment + (1 |", rand, ")")), b, REML = TRUE)
  vc <- as.data.frame(VarCorr(m))
  # emmeans applies Kenward-Roger to lmer fits by default, which inflates the SEs; the GLS (plug-in)
  # SEs of recover_interblock() correspond to lmer.df = "asymptotic"
  em <- summary(emmeans(m, "Treatment", lmer.df = "asymptotic"))
  em <- em[match(rc$means$level, em$Treatment), ]
  chk(paste(label, "REML variance components vs lme4"), near(c(rc$sigma2_block, rc$sigma2), vc$vcov, 1e-4),
      sprintf("(%.4f, %.4f)", rc$sigma2_block, rc$sigma2))
  chk(paste(label, "combined treatment means and SEs vs lme4 + emmeans"), near(rc$means$mean, em$emmean, 1e-4) && near(rc$means$se, em$SE, 1e-3))
  rc
}
rc <- rec_chk("alpha t=24 k=4 r=3:", sim(make_design("alpha", trt = 24, k = 4, r = 3, seed = 7), 1.2), "Block:IncBlock", block = c("Block", "IncBlock"))
print(rc$sed)
rec_chk("lattice 5 x 5 r=2:", sim(make_design("lattice", trt = 25, k = 5, r = 2, seed = 7), 1), "Block:IncBlock", block = c("Block", "IncBlock"))
rec_chk("BIBD t=7 k=3:", sim(make_design("bib", trt = 7, k = 3, seed = 7), 1), "Block", block = "Block")
b <- sim(make_design("augmented", checks = c("C1", "C2", "C3"), entries = paste0("E", 1:15), r = 4, seed = 7), 1)
rec_chk("augmented:", b, "Block", block = "Block", design = "augmented")
cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
