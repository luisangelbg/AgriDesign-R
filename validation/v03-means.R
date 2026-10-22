# Stage 3 validation: mean separation, Dunnett, simple effects, trends and
# contrasts against base R, agricolae, multcomp, emmeans, lmerTest, PMCMRplus
# and the figures of the web application (user manual, chapter 6).
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages({ library(agricolae); library(multcomp); library(emmeans); library(lmerTest); library(PMCMRplus) })
ok <- 0; bad <- 0
chk <- function(what, got, ref, tol = 1e-6) {
  got <- as.numeric(got); ref <- as.numeric(ref)
  d <- max(abs(got - ref) / pmax(1, abs(ref)))
  pass <- length(got) == length(ref) && is.finite(d) && d < tol
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-60s got %-24s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(signif(got, 6), collapse = " "), paste(signif(ref, 6), collapse = " ")))
}
same <- function(what, got, ref) {
  pass <- identical(as.character(got), as.character(ref))
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-60s got %-24s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(got, collapse = " "), paste(ref, collapse = " ")))
}
# p of the pair (a, b) from a lower-triangular matrix of p-values
pm_get <- function(M, a, b) mapply(function(x, y) { v <- M[x, y]; if (is.na(v)) M[y, x] else v }, a, b)
# letters of agricolae re-ordered as our means
agl <- function(ag, labels) as.character(ag$groups[as.character(labels), "groups"])

# ================= maize (RCBD) =================
mz <- maize_nitrogen; mz$N <- factor(mz$Nitrogen_kg_ha); mz$B <- factor(mz$Block)
av <- aov(Yield_t_ha ~ B + N, mz)
f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
tk <- compare_means(f, "Nitrogen_kg_ha")
th <- TukeyHSD(av, "N")$N
chk("maize Tukey: p vs TukeyHSD", tk$pairs$p, th[paste(tk$pairs$level_j, tk$pairs$level_i, sep = "-"), "p adj"], 1e-6)
chk("maize Tukey: HSD vs app 0.283", round(tk$msd, 3), 0.283)
chk("maize Tukey: 120 vs 180 p vs app 0.0393", round(tk$pairs$p[tk$pairs$level_i == "120" & tk$pairs$level_j == "180"], 4), 0.0393)
same("maize Tukey: letters d c b a (0, 60, 120, 180)", tk$means$group, c("d", "c", "b", "a"))
ag <- HSD.test(av, "N"); chk("maize Tukey: HSD vs agricolae", tk$msd, ag$statistics$MSD)
l <- compare_means(f, "Nitrogen_kg_ha", method = "lsd"); ag <- LSD.test(av, "N")
chk("maize LSD vs agricolae LSD.test", l$msd, ag$statistics$LSD)
d <- compare_means(f, "Nitrogen_kg_ha", method = "duncan"); ag <- duncan.test(av, "N")
chk("maize Duncan: critical ranges vs agricolae", d$msd, ag$duncan$CriticalRange)
same("maize Duncan: letters vs agricolae", d$means$group, agl(ag, d$means$label))
s <- compare_means(f, "Nitrogen_kg_ha", method = "snk"); ag <- SNK.test(av, "N")
chk("maize SNK: critical ranges vs agricolae", s$msd, ag$snk$CriticalRange)
same("maize SNK: letters vs agricolae", s$means$group, agl(ag, s$means$label))
r <- compare_means(f, "Nitrogen_kg_ha", method = "regwq"); ag <- REGW.test(av, "N")
chk("maize REGWQ: critical ranges vs agricolae", r$msd, ag$regw$CriticalRange)
same("maize REGWQ: letters vs agricolae", r$means$group, agl(ag, r$means$label))
sc <- compare_means(f, "Nitrogen_kg_ha", method = "scheffe"); ag <- scheffe.test(av, "N")
chk("maize Scheffe: critical difference vs agricolae", sc$msd, ag$statistics$CriticalDifference)
same("maize Scheffe: letters vs agricolae", sc$means$group, agl(ag, sc$means$label))
# trend (manual: SS linear 19.7707, quadratic 1.7756; y = 4.207 + 0.03323 x - 0.00009253 x^2, max at 179.5)
tr <- trend(f, "Nitrogen_kg_ha")
chk("maize trend: SS linear, quadratic vs app", round(tr$components$ss[1:2], 4), c(19.7707, 1.7756))
avp <- aov(Yield_t_ha ~ B + N, mz, contrasts = list(N = contr.poly(4)))
sp <- summary(avp, split = list(N = list(L = 1, Q = 2, C = 3)))[[1]]
chk("maize trend: SS L, Q, C vs aov split contr.poly", tr$components$ss, sp[c(" N: L", " N: Q", " N: C"), "Sum Sq"] |> unlist() |>
      (\(v) if (anyNA(v)) sp[grep("N: ", rownames(sp)), "Sum Sq"] else v)())
chk("maize trend: quadratic coefficients vs app", signif(tr$coef, 4), c(4.207, 0.03323, -0.00009253), 1e-3)
chk("maize trend: maximum at 179.5", round(tr$optimum[["x"]], 1), 179.5)
chk("maize trend: polynomial vs lm(poly raw)", tr$coef, coef(lm(Yield_t_ha ~ poly(Nitrogen_kg_ha, 2, raw = TRUE), mz)))
# contrast (manual: planned 120 vs 180 p = 0.0095)
ct <- test_contrasts(f, "Nitrogen_kg_ha", rbind("120 vs 180" = c(0, 0, 1, -1), "0 vs rest" = c(3, -1, -1, -1)))
chk("maize contrast 120 vs 180: p vs app 0.0095", round(ct$p[1], 4), 0.0095)
em <- summary(contrast(emmeans(av, "N"), list(a = c(0, 0, 1, -1), b = c(3, -1, -1, -1))))
chk("maize contrasts: estimate, SE, p vs emmeans", c(ct$estimate, ct$se, ct$p), c(em$estimate, em$SE, em$p.value))

# ================= beans (CRD) =================
bv <- bean_varieties
f <- design_anova(bv, "Yield_kg_plot", "Variety")
lv <- levels(f$data$Variety)
bv$V <- factor(bv$Variety, levels = lv)
av <- aov(Yield_kg_plot ~ V, bv)
tk <- compare_means(f, "Variety")
chk("beans Tukey: HSD vs app 0.224", round(tk$msd, 3), 0.224)
bo <- compare_means(f, "Variety", method = "bonferroni")
chk("beans Bonferroni: MSD vs app 0.239", round(bo$msd, 3), 0.239)
pt <- pairwise.t.test(bv$Yield_kg_plot, bv$V, p.adjust.method = "bonferroni")$p.value
chk("beans Bonferroni: p vs pairwise.t.test", bo$pairs$p_adj, pm_get(pt, bo$pairs$level_j, bo$pairs$level_i))
ho <- compare_means(f, "Variety", method = "holm")
pt <- pairwise.t.test(bv$Yield_kg_plot, bv$V, p.adjust.method = "holm")$p.value
chk("beans Holm: p vs pairwise.t.test", ho$pairs$p_adj, pm_get(pt, ho$pairs$level_j, ho$pairs$level_i))
si <- compare_means(f, "Variety", method = "sidak")
pt <- pairwise.t.test(bv$Yield_kg_plot, bv$V, p.adjust.method = "none")$p.value
chk("beans Sidak: p = 1 - (1 - p)^m", si$pairs$p_adj, 1 - (1 - pm_get(pt, si$pairs$level_j, si$pairs$level_i))^10)
# Dunnett (manual: critical difference 0.198, Bayo vs Negro Jamapa p = 0.0318)
dn <- compare_means(f, "Variety", method = "dunnett", control = "Negro Jamapa")
chk("beans Dunnett: critical difference vs app 0.198", round(dn$msd, 3), 0.198)
chk("beans Dunnett: Bayo vs Negro Jamapa p vs app", round(dn$pairs$p[dn$pairs$level_i == "Bayo Victoria" | grepl("^Bayo", dn$pairs$level_i)], 4), 0.0318)
# Dunnett reference: mvtnorm at high precision (multcomp's default tolerance, 1e-3, is too coarse:
# it gives a critical t of 2.7251 where the exact value is 2.7273). Very small p-values need
# many points: with 2e6 points mvtnorm gives 1.96e-6 for t = 9.04 (unbalanced case), with 2e7
# 3.89e-6 (error 4e-8); ours is 3.93e-6, inside the exact bounds [p single, 4 p single].
mvt_p <- function(t, n, n0, df) {
  lam <- sqrt(n / (n + n0)); R <- outer(lam, lam); diag(R) <- 1
  set.seed(1)
  vapply(abs(t), function(tt) 1 - mvtnorm::pmvt(lower = rep(-tt, length(n)), upper = rep(tt, length(n)), corr = R, df = df,
                                                algorithm = mvtnorm::GenzBretz(maxpts = 2e7, abseps = 1e-12))[1], 0)
}
nj <- f$data$Variety |> table()
ref <- mvt_p(dn$pairs$diff / dn$pairs$se, as.numeric(nj[dn$pairs$level_i]), nj[["Negro Jamapa"]], f$residual[["df"]])
chk("beans Dunnett: p vs mvtnorm (high precision)", dn$pairs$p, ref, 1e-6)
tc <- dn$msd / dn$pairs$se[1]
chk("beans Dunnett: p at our critical t = 0.05 (mvtnorm)", mvt_p(tc, rep(4, 4), 4, f$residual[["df"]]), 0.05, 1e-6)
# unequal replication: Tukey-Kramer and exact Dunnett correlations
bu <- bv[-c(2, 7, 9), ]
fu <- design_anova(bu, "Yield_kg_plot", "Variety")
bu$V <- factor(bu$Variety, levels = levels(fu$data$Variety))
tku <- compare_means(fu, "Variety")
thu <- TukeyHSD(aov(Yield_kg_plot ~ V, bu), "V")$V
chk("unbalanced Tukey-Kramer: p vs TukeyHSD", tku$pairs$p, thu[paste(tku$pairs$level_j, tku$pairs$level_i, sep = "-"), "p adj"], 1e-6)
dnu <- compare_means(fu, "Variety", method = "dunnett", control = "Negro Jamapa")
nu <- table(fu$data$Variety)
ref <- mvt_p(dnu$pairs$diff / dnu$pairs$se, as.numeric(nu[dnu$pairs$level_i]), nu[["Negro Jamapa"]], fu$residual[["df"]])
chk("unbalanced Dunnett: p vs mvtnorm (exact correlations)", dnu$pairs$p, ref, 1e-6)
gh <- compare_means(f, "Variety", method = "games_howell")
ghp <- gamesHowellTest(bv$Yield_kg_plot, bv$V)$p.value
chk("beans Games-Howell: p vs PMCMRplus", gh$pairs$p, pm_get(ghp, gh$pairs$level_j, gh$pairs$level_i), 1e-6)
d <- compare_means(f, "Variety", method = "duncan"); ag <- duncan.test(av, "V")
same("beans Duncan: letters vs agricolae", d$means$group, agl(ag, d$means$label))

# ================= ANCOVA: Tukey on adjusted means vs emmeans =================
cv <- maize_nitrogen; cv$Stand <- c(58, 61, 63, 57, 60, 62, 59, 64, 55, 61, 60, 62, 63, 58, 59, 60)
f <- design_anova(cv, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block", covariate = "Stand")
tk <- compare_means(f, "Nitrogen_kg_ha")
cv$N <- factor(cv$Nitrogen_kg_ha); cv$B <- factor(cv$Block)
pe <- summary(pairs(emmeans(lm(Yield_t_ha ~ Stand + B + N, cv), "N"), reverse = FALSE))
chk("ANCOVA Tukey: SE of differences vs emmeans", tk$pairs$se, pe$SE)
chk("ANCOVA Tukey: p vs emmeans", tk$pairs$p, pe$p.value, 1e-5)

# ================= split-plot: strata and simple effects =================
sp <- split_irrigation; sp$Block <- factor(sp$Block)
f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd")
mi <- compare_means(f, "Irrigation", method = "lsd")
spa <- with(sp, sp.plot(Block, Irrigation, Variety, Yield_t_ha))
ag <- with(sp, LSD.test(Yield_t_ha, Irrigation, spa$gl.a, spa$Ea))
chk("split: LSD of Irrigation (Error a) vs agricolae", mi$msd, ag$statistics$LSD)
ag <- with(sp, LSD.test(Yield_t_ha, Variety, spa$gl.b, spa$Eb))
chk("split: LSD of Variety (Error b) vs agricolae", compare_means(f, "Variety", method = "lsd")$msd, ag$statistics$LSD)
mm <- lmer(Yield_t_ha ~ Irrigation * Variety + (1 | Block) + (1 | Block:Irrigation), sp)
se1 <- compare_means(f, "Variety", by = "Irrigation", method = "lsd")
e1 <- summary(pairs(emmeans(mm, ~ Variety | Irrigation, lmer.df = "satterthwaite"), adjust = "none"))
chk("split: Variety within Irrigation, SE diff vs lmer", se1$pairs$se, e1$SE, 1e-5)
chk("split: Variety within Irrigation, df vs lmer", rep(se1$error$df, nrow(e1)), e1$df, 1e-4)
se2 <- compare_means(f, "Irrigation", by = "Variety", method = "lsd")
e2 <- summary(pairs(emmeans(mm, ~ Irrigation | Variety, lmer.df = "satterthwaite"), adjust = "none"))
chk("split: Irrigation within Variety, SE diff vs lmer", se2$pairs$se, e2$SE, 1e-5)
chk("split: Irrigation within Variety, Satterthwaite df vs lmer", rep(se2$error$df, nrow(e2)), e2$df, 1e-4)
chk("split: p of Irrigation within Variety vs lmer", se2$pairs$p, e2$p.value, 1e-5)
# strip-plot: simple effects with two pooled errors
f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "strip_rcbd")
ms <- lmer(Yield_t_ha ~ Irrigation * Variety + (1 | Block) + (1 | Block:Irrigation) + (1 | Block:Variety), sp)
s1 <- compare_means(f, "Variety", by = "Irrigation", method = "lsd")
e1 <- summary(pairs(emmeans(ms, ~ Variety | Irrigation, lmer.df = "satterthwaite"), adjust = "none"))
chk("strip: Variety within Irrigation, SE and df vs lmer", c(s1$pairs$se, s1$error$df), c(e1$SE, e1$df[1]), 1e-4)
s2 <- compare_means(f, "Irrigation", by = "Variety", method = "lsd")
e2 <- summary(pairs(emmeans(ms, ~ Irrigation | Variety, lmer.df = "satterthwaite"), adjust = "none"))
chk("strip: Irrigation within Variety, SE and df vs lmer", c(s2$pairs$se, s2$error$df), c(e2$SE, e2$df[1]), 1e-4)

# ================= factorial: all cells together, back-transformed means =================
f <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block", transform = "log")
cc <- compare_means(f, c("Variety", "Fertilizer"))
tm <- tomato_factorial; tm$B <- factor(tm$Block)
th <- TukeyHSD(aov(log(Yield_kg_plant) ~ B + Variety * Fertilizer, tm), "Variety:Fertilizer")[[1]]
key <- function(l) sub(" x ", ":", l)
thp <- function(a, b) { n1 <- paste(key(a), key(b), sep = "-"); n2 <- paste(key(b), key(a), sep = "-")
  ifelse(n1 %in% rownames(th), th[match(n1, rownames(th)), "p adj"], th[match(n2, rownames(th)), "p adj"]) }
chk("tomato cells (log): Tukey p vs TukeyHSD", cc$pairs$p, thp(cc$pairs$level_i, cc$pairs$level_j), 1e-6)
chk("tomato: back-transformed means = exp(LS mean)", cc$means$back_transformed, exp(cc$means$mean))

cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
