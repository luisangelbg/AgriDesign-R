# Stage 7c validation: Waller-Duncan k-ratio t test, Durbin, van der Waerden and
# Mood's median tests, AUDPC / AUDPS and Smith's index of soil heterogeneity,
# against agricolae, PMCMRplus, base R, SAS and independent computations.
# agricolae's data sets are used here only (never inside the package).
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages({ library(agricolae); library(PMCMRplus) })
# agricolae also exports audpc() and audps(): ours are called as AgriDesign:::audpc() below
# S3 methods of the new files, until roxygen writes them into NAMESPACE
for (m in c("print.ad_durbin", "print.ad_waerden", "print.ad_median_test", "print.ad_smith", "plot.ad_smith"))
  registerS3method(sub("\\..*", "", m), sub("^[a-z]+\\.", "", m), get(m, asNamespace("AgriDesign")))
ok <- 0; bad <- 0
chk <- function(what, got, ref, tol = 1e-6) {
  got <- as.numeric(got); ref <- as.numeric(ref)
  d <- max(abs(got - ref) / pmax(1, abs(ref)))
  pass <- length(got) == length(ref) && is.finite(d) && d < tol
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-62s got %-28s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(signif(got, 6), collapse = " "), paste(signif(ref, 6), collapse = " ")))
}
same <- function(what, got, ref) {
  pass <- identical(as.character(got), as.character(ref))
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-62s got %-28s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(got, collapse = " "), paste(ref, collapse = " ")))
}
note <- function(...) cat("NOTE", sprintf(...), "\n")
# p-value of the pair (a, b) in a PMCMRplus lower-triangular matrix
pm_get <- function(P, a, b) mapply(function(i, j) { v <- P[i, j]; if (is.na(v)) P[j, i] else v }, a, b)
# letters of agricolae's groups, aligned with our levels
ag_letters <- function(gr, lv) as.character(gr$groups[match(lv, rownames(gr))])

# ======================= 1. Waller-Duncan =======================
cat("\n--- Waller-Duncan critical t ---\n")
# (a) SAS/STAT ANOVA documentation, Example "Multiple comparisons" (red clover,
#     6 strains x 5 plots): K = 100, error df 24, F = 14.37, critical t 1.91873,
#     MSE 11.78867, minimum significant difference 4.1665.
tb <- waller_critical(100, 5, 24, 14.37)
chk("SAS red clover: t_B(100, q=5, f=24, F=14.37)", round(tb, 5), 1.91873, 1e-9)
chk("SAS red clover: MSD = t_B sqrt(2 MSE / 5) (printed F, MSE rounded)", tb * sqrt(2 * 11.78867 / 5), 4.1665, 2e-5)

# (b) independent evaluation: the posterior expected losses as a double integral
#     over the two chi-square variates (no closed-form reduction), solved for t.
g0 <- function(z) z * pnorm(z) + dnorm(z)
bf_E <- function(t, q, f, Fv, sgn) {
  # both chi-square variates on the probability scale (smooth integrands)
  outer_f <- function(pa) vapply(pa, function(p) {
    ai <- qchisq(p, f); u <- f / ai; pmax_b <- pchisq(q * ai * Fv / f, q)
    inner <- function(pb) { b <- qchisq(pb, q); beta <- pmax(0, 1 - f * b / (q * ai * Fv)); sqrt(u * beta) * g0(sgn * t * sqrt(beta / u)) }
    integrate(inner, 0, pmax_b, rel.tol = 1e-8, subdivisions = 500, stop.on.error = FALSE)$value
  }, 0)
  integrate(outer_f, 0, 1, rel.tol = 1e-8, subdivisions = 500, stop.on.error = FALSE)$value
}
bf_t <- function(K, q, f, Fv) {
  near <- waller_critical(K, q, f, Fv)  # only brackets the root; the root itself comes from the double integral
  uniroot(function(t) log(bf_E(t, q, f, Fv, 1) / bf_E(t, q, f, Fv, -1)) - log(K), near + c(-0.05, 0.05), tol = 1e-7, extendInt = "yes")$root
}
grid <- rbind(c(1, 20, 5), c(1, 5, 2), c(2, 10, 3), c(4, 20, 10), c(10, 20, 100), c(3, 8, 1.5))
for (i in seq_len(nrow(grid))) {
  q <- grid[i, 1]; f <- grid[i, 2]; Fv <- grid[i, 3]
  chk(sprintf("t_B vs double integral, q=%g f=%g F=%g", q, f, Fv), waller_critical(100, q, f, Fv), bf_t(100, q, f, Fv), 1e-6)
}
chk("t_B K=500 vs double integral, q=4 f=12 F=6", waller_critical(500, 4, 12, 6), bf_t(500, 4, 12, 6), 1e-6)
chk("t_B K=50 vs double integral, q=6 f=30 F=4", waller_critical(50, 6, 30, 4), bf_t(50, 6, 30, 4), 1e-6)

# (c) limits: F -> Inf gives the t-based rule (independent of q); f, F -> Inf the normal rule.
lim_t <- function(K, f) uniroot(function(t) {
  e <- function(s) integrate(function(a) sqrt(f / a) * g0(s * t * sqrt(a / f)) * dchisq(a, f), 0, Inf, rel.tol = 1e-12)$value
  log(e(1) / e(-1)) - log(K) }, c(0.5, 6), tol = 1e-10)$root
chk("F -> Inf (1e8), q = 1, 4, 10, f = 20 vs t-based limit", waller_critical(100, c(1, 4, 10), 20, 1e8), rep(lim_t(100, 20), 3), 1e-5)
zl <- uniroot(function(z) g0(z) / g0(-z) - 100, c(0.5, 5), tol = 1e-12)$root
chk("f, F -> Inf vs normal limit g(z)/g(-z) = 100 (1.7208)", waller_critical(100, 3, 1e7, 1e8), zl, 1e-4)

# (d) agricolae:::waller (Waller-Kemp algorithm, printed to 3 decimals)
cat("\n--- t_B vs agricolae:::waller ---\n")
cmp <- expand.grid(F = c(1.5, 2, 3, 5, 10, 100), q = c(1, 4, 10), f = c(5, 20, 60))
cmp$mine <- waller_critical(100, cmp$q, cmp$f, cmp$F)
cmp$agricolae <- mapply(function(q, f, Fv) agricolae:::waller(100, q, f, Fv), cmp$q, cmp$f, cmp$F)
cmp$diff <- cmp$mine - cmp$agricolae
agree <- abs(cmp$diff) <= 0.0011
chk(sprintf("t_B agrees with agricolae to 3 decimals in %d of %d cells with q >= 2, F <= 10",
            sum(agree & cmp$q > 1 & cmp$F <= 10), sum(cmp$q > 1 & cmp$F <= 10)),
    sum(agree & cmp$q > 1 & cmp$F <= 10), sum(cmp$q > 1 & cmp$F <= 10), 1e-12)
dis <- cmp[!agree, ]
note("agricolae differs (> 0.001) in %d cells (q = 1, or F = 100):", nrow(dis))
print(dis, row.names = FALSE, digits = 5)
sel <- dis[c(1, 6, 8, 9, 16, 24), ]
sel$double_integral <- mapply(function(q, f, Fv) bf_t(100, q, f, Fv), sel$q, sel$f, sel$F)
chk("6 of these cells: AgriDesign = independent double integral", sel$mine, sel$double_integral, 1e-6)
note("agricolae is off by up to %.3f there: its Waller-Kemp series is inexact for q = 1 and for large F", max(abs(sel$agricolae - sel$double_integral)))

# (e) waller.test on agricolae's sweetpotato (CRD, 4 viruses x 3): q = 3, f = 8, F = 17.34
cat("\n--- waller_duncan() vs waller.test ---\n")
data(sweetpotato)
w_ag <- waller.test(aov(yield ~ virus, sweetpotato), "virus", K = 100)
fit <- design_anova(sweetpotato, "yield", "virus")
wd <- waller_duncan(fit, "virus")
chk("sweetpotato: F", fit$table$F[1], w_ag$statistics$F.Value, 1e-6)
chk("sweetpotato: t_B (agricolae rounds to 3 decimals)", round(wd$t_critical, 3), w_ag$statistics$Waller, 1e-12)
chk("sweetpotato: t_B vs double integral", wd$t_critical, bf_t(100, 3, 8, fit$table$F[1]), 1e-6)
chk("sweetpotato: MSD (agricolae uses the rounded t)", wd$msd / wd$t_critical * round(wd$t_critical, 3),
    w_ag$statistics$CriticalDifference, 1e-8)
same("sweetpotato: letters", wd$means$group, ag_letters(w_ag$groups, wd$means$label))
wd5 <- waller_duncan(fit, "virus", alpha = 0.01)
same("alpha = 0.01 sets K = 500", wd5$K, 500)
# maize (RCBD): F from the ANOVA of the package, f = 9
fm <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
wm <- waller_duncan(fm, "Nitrogen_kg_ha")
am <- aov(Yield_t_ha ~ factor(Block) + factor(Nitrogen_kg_ha), maize_nitrogen)
w_agm <- waller.test(am, "factor(Nitrogen_kg_ha)", K = 100)
chk("maize RCBD (F = 437): t_B vs double integral", wm$t_critical, bf_t(100, 3, 9, fm$table$F[2]), 1e-6)
note("maize RCBD: agricolae t_B = %.3f, AgriDesign %.4f (large F, see above)", w_agm$statistics$Waller, wm$t_critical)
same("maize RCBD: letters vs agricolae", wm$means$group, ag_letters(w_agm$groups, wm$means$label))
# split-plot: the main-plot factor uses Error a for F, f and the SE
fs <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd")
ws <- waller_duncan(fs, "Variety")
ra <- fs$table[fs$table$term == "Variety", ]
chk("split-plot sub-plot factor: t_B from its own F and df", ws$t_critical, waller_critical(100, ra$df, ra$df_error, ra$F))

# ======================= 2. Durbin =======================
cat("\n--- Durbin ---\n")
bs <- bib_sorghum
du <- durbin_test(bs, "Yield_t_ha", "Hybrid", "Block", p_adjust = "none")
pd <- durbinTest(bs$Yield_t_ha, bs$Hybrid, bs$Block)
chk("bib_sorghum: chi-square vs PMCMRplus::durbinTest", du$test$statistic, pd$statistic, 1e-10)
chk("bib_sorghum: p vs PMCMRplus::durbinTest", du$test$p, pd$p.value, 1e-10)
pa <- durbinAllPairsTest(bs$Yield_t_ha, bs$Hybrid, bs$Block, p.adjust.method = "none")
chk("bib_sorghum: pairwise p vs durbinAllPairsTest", du$pairs$p, pm_get(pa$p.value, du$pairs$level_j, du$pairs$level_i), 1e-8)
duh <- durbin_test(bs, "Yield_t_ha", "Hybrid", "Block")
pah <- durbinAllPairsTest(bs$Yield_t_ha, bs$Hybrid, bs$Block, p.adjust.method = "holm")
chk("bib_sorghum: Holm p vs durbinAllPairsTest(holm)", duh$pairs$p_adj, pm_get(pah$p.value, duh$pairs$level_j, duh$pairs$level_i), 1e-8)
judge <- bs$Block; trt <- bs$Hybrid; evaluation <- bs$Yield_t_ha
da <- durbin.test(judge, trt, evaluation, alpha = 0.05, group = TRUE)
chk("bib_sorghum: chi-square, p vs agricolae", c(du$test$statistic, du$test$p), c(da$statistics$chisq.value, da$statistics$p.value), 1e-8)
chk("bib_sorghum: critical t of the LSD vs agricolae", qt(0.975, du$df_error), da$statistics$t.value, 1e-8)
chk("bib_sorghum: LSD of rank sums vs agricolae", du$lsd_ranks, da$statistics$LSD, 1e-8)
same("bib_sorghum: letters vs agricolae", du$groups$group, ag_letters(da$groups, du$groups$level))
# Conover (1999, p. 391) ice-cream example, as in the help of agricolae::durbin.test
runs <- gl(7, 3)
ic <- c(1, 2, 4, 2, 3, 5, 3, 4, 6, 4, 5, 7, 1, 5, 6, 2, 6, 7, 1, 3, 7)
val <- c(2, 3, 1, 3, 1, 2, 2, 1, 3, 1, 2, 3, 3, 1, 2, 3, 1, 2, 3, 1, 2)
dc <- durbin_test(data.frame(r = runs, t = ic, y = val), "y", "t", "r", p_adjust = "none")
judge <- runs; trt <- ic; evaluation <- val
dca <- durbin.test(judge, trt, evaluation, group = TRUE)
chk("Conover ice cream: chi-square vs agricolae", dc$test$statistic, dca$statistics$chisq.value, 1e-10)
chk("Conover ice cream: chi-square vs PMCMRplus::durbinTest", dc$test$statistic, durbinTest(val, factor(ic), runs)$statistic, 1e-10)
chk("Conover ice cream: LSD of rank sums vs agricolae", dc$lsd_ranks, dca$statistics$LSD, 1e-8)
# complete blocks: Durbin = Friedman
dsf <- durbin_test(design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block"))
chk("complete blocks: Durbin = friedman.test", dsf$test$statistic,
    friedman.test(Severity_score ~ Fungicide | Block, disease_scores)$statistic, 1e-10)

# ======================= 3. van der Waerden =======================
cat("\n--- van der Waerden ---\n")
data(corn)
cn <- data.frame(y = corn$observation, m = factor(corn$method))
vw <- waerden_test(cn, "y", "m", p_adjust = "none")
pv <- vanWaerdenTest(cn$y, cn$m)
chk("corn: statistic vs PMCMRplus::vanWaerdenTest", vw$test$statistic, pv$statistic, 1e-10)
pva <- vanWaerdenAllPairsTest(cn$y, cn$m, p.adjust.method = "none")
chk("corn: pairwise p vs vanWaerdenAllPairsTest", vw$pairs$p, pm_get(pva$p.value, vw$pairs$level_j, vw$pairs$level_i), 1e-8)
wa <- waerden.test(cn$y, cn$m, group = FALSE)
note("corn: agricolae chi-square %.4f vs %.4f here (and PMCMRplus): agricolae rounds R/(N + 1) to 3 decimals before qnorm()", wa$statistics$Chisq, vw$test$statistic)
Nc <- nrow(cn); Ar <- qnorm(round(rank(cn$y) / (Nc + 1), 3)); Abr <- tapply(Ar, cn$m, mean)
T1r <- sum(table(cn$m) * Abr^2) / (sum(Ar^2) / (Nc - 1))
chk("corn: agricolae chi-square = same formula with rounded probabilities", T1r, wa$statistics$Chisq, 1e-10)
chk("corn: agricolae normal scores = rounded-probability scores", Abr, wa$means$normalScore, 1e-10)
chk("corn: mean normal scores = exact qnorm(R / (N + 1)) means", vw$groups$mean_score[order(vw$groups$level)], tapply(qnorm(rank(cn$y) / (Nc + 1)), cn$m, mean), 1e-12)
chk("corn: pairwise p vs agricolae to 3 decimals (rounding of scores)", round(vw$pairs$p, 3), round(wa$comparison$pvalue, 3), 1e-12)
wg <- waerden.test(cn$y, cn$m, group = TRUE)
same("corn: letters vs agricolae", vw$groups$group, ag_letters(wg$groups, vw$groups$level))
vb <- waerden_test(design_anova(bean_varieties, "Yield_kg_plot", "Variety"))
chk("beans (fit): statistic vs vanWaerdenTest", vb$test$statistic, vanWaerdenTest(bean_varieties$Yield_kg_plot, factor(bean_varieties$Variety))$statistic, 1e-10)

# ======================= 4. Mood's median test =======================
cat("\n--- Median test ---\n")
mt <- median_test(cn, "y", "m", p_adjust = "none")
tab <- table(cn$y > median(cn$y), cn$m)
cs <- suppressWarnings(chisq.test(tab, correct = FALSE))
chk("corn: overall chi-square, p vs chisq.test on the 2 x k table", c(mt$test$statistic, mt$test$p), c(cs$statistic, cs$p.value), 1e-10)
chk("corn: overall chi-square vs agricolae", mt$test$statistic, (ma <- Median.test(cn$y, cn$m, console = FALSE, group = FALSE))$statistics$Chisq, 1e-10)
pp <- vapply(seq_len(nrow(mt$pairs)), function(q) {
  s <- cn$m %in% c(mt$pairs$level_i[q], mt$pairs$level_j[q]); yy <- cn$y[s]
  suppressWarnings(chisq.test(table(yy > median(yy), droplevels(cn$m[s])), correct = FALSE)$p.value) }, 0)
chk("corn: pairwise p vs chisq.test on each 2 x 2 (pair median)", mt$pairs$p, pp, 1e-10)
ag_p <- ma$comparison$pvalue
chk("corn: pairwise p vs agricolae (rounded to 4 decimals)", round(mt$pairs$p, 4), ag_p, 1e-12)
mtc <- median_test(cn, "y", "m", p_adjust = "none", correct = TRUE)
ppc <- vapply(seq_len(nrow(mt$pairs)), function(q) {
  s <- cn$m %in% c(mt$pairs$level_i[q], mt$pairs$level_j[q]); yy <- cn$y[s]
  suppressWarnings(chisq.test(table(yy > median(yy), droplevels(cn$m[s])), correct = TRUE)$p.value) }, 0)
chk("corn: Yates-corrected pairwise p vs chisq.test(correct = TRUE)", mtc$pairs$p, ppc, 1e-10)
mag_c <- Median.test(cn$y, cn$m, correct = TRUE, console = FALSE, group = FALSE)$comparison$pvalue
note("agricolae Median.test(correct = TRUE) p = %s equal the UNcorrected ones: chisq.test(B, correct, ...) passes 'correct' as 'y'",
     paste(head(mag_c, 3), collapse = " "))
chk("agricolae's 'corrected' pairwise p = uncorrected chi-square (its bug)", mag_c, round(pp, 4), 1e-12)
mte <- median_test(cn, "y", "m", p_adjust = "none", exact = TRUE)
pe <- vapply(seq_len(nrow(mt$pairs)), function(q) {
  s <- cn$m %in% c(mt$pairs$level_i[q], mt$pairs$level_j[q]); yy <- cn$y[s]
  fisher.test(table(yy > median(yy), droplevels(cn$m[s])))$p.value }, 0)
chk("corn: exact pairwise p vs fisher.test", mte$pairs$p, pe, 1e-10)
mg <- median_test(cn, "y", "m", p_adjust = "none", pair_median = "grand")
pg <- vapply(seq_len(nrow(mt$pairs)), function(q) {
  s <- cn$m %in% c(mt$pairs$level_i[q], mt$pairs$level_j[q])
  suppressWarnings(chisq.test(table(factor(cn$y[s] > median(cn$y), c(FALSE, TRUE)), droplevels(cn$m[s])), correct = FALSE)$p.value) }, 0)
chk("corn: grand-median pairwise p vs chisq.test", mg$pairs$p, pg, 1e-10)
mg2 <- Median.test(cn$y, cn$m, console = FALSE, group = TRUE)
same("corn: letters vs agricolae", mt$groups$group, ag_letters(mg2$groups, mt$groups$level))
note("mood.test() is a test of scale, not the median test, so it is not a reference here")

# ======================= 5. AUDPC / AUDPS =======================
cat("\n--- AUDPC and AUDPS ---\n")
dates <- c(14, 21, 28, 35, 42)
ev <- data.frame(E1 = c(40, 80, 10), E2 = c(80, 90, 30), E3 = c(90, 98, 32), E4 = c(95, 99, 40), E5 = c(100, 100, 63))
chk("AUDPC (3 units) vs agricolae::audpc", AgriDesign:::audpc(ev, dates), agricolae::audpc(ev, dates), 1e-12)
chk("relative AUDPC vs agricolae (both / (100 D))", AgriDesign:::audpc(ev, dates, "relative"), agricolae::audpc(ev, dates, "relative"), 1e-12)
chk("AUDPS (data frame) vs agricolae::audps", AgriDesign:::audps(ev, dates), unlist(agricolae::audps(ev, dates)), 1e-12)
evm <- as.matrix(ev)
ag_m <- agricolae::audps(evm, dates)
note("agricolae::audps on a MATRIX uses evaluation[1] and evaluation[n] (single cells): %s vs %s",
     paste(round(ag_m, 1), collapse = " "), paste(round(AgriDesign:::audps(evm, dates), 1), collapse = " "))
chk("AUDPS of each matrix row vs agricolae one row at a time", AgriDesign:::audps(evm, dates),
    vapply(1:3, function(i) agricolae::audps(evm[i, ], dates), 0), 1e-12)
n <- length(dates)
chk("relative AUDPS x n/(n-1) = agricolae's (it divides by 100 D only)", AgriDesign:::audps(ev, dates, "relative") * n / (n - 1),
    unlist(agricolae::audps(ev, dates, "relative")), 1e-12)
note("agricolae's relative AUDPS reaches n/(n-1) = %.2f when y = 100 always; ours reaches 1 (Simko & Piepho 2012)", n / (n - 1))
chk("relative AUDPS of y = ymax always is 1", AgriDesign:::audps(rep(100, 5), dates, "relative"), 1)
chk("relative AUDPC of y = ymax always is 1", AgriDesign:::audpc(rep(100, 5), dates, "relative"), 1)
# independent identities
dd <- c(0, 3, 10, 12, 20, 31); yy <- c(0, 2, 7, 15, 33, 60)
lin <- approxfun(dd, yy)
chk("AUDPC = integral of the linear interpolation (integrate by segment)", AgriDesign:::audpc(yy, dd), sum(mapply(function(a, b) integrate(lin, a, b, rel.tol = 1e-12)$value, dd[-6], dd[-1])), 1e-10)
chk("AUDPS with equal intervals = sum(y) * interval", AgriDesign:::audps(ev, dates), rowSums(ev) * 7, 1e-12)
chk("standardised AUDPC of constant y is y", AgriDesign:::audpc(rep(37, 4), c(1, 5, 6, 20), "standardized"), 37)
chk("ymax = 1 (proportions) gives the same relative areas", AgriDesign:::audps(ev / 100, dates, "relative", ymax = 1), AgriDesign:::audps(ev, dates, "relative"), 1e-12)
chk("Date values give the same areas as day numbers", AgriDesign:::audpc(yy, as.Date("2026-05-01") + dd), AgriDesign:::audpc(yy, dd), 1e-12)
# per plot from the long data
au <- audpc_by_unit(blight_progress, "Severity_pct", "Day", c("Block", "Treatment"))
w <- reshape(blight_progress, idvar = c("Block", "Treatment"), timevar = "Day", direction = "wide")
w <- w[match(paste(au$Block, au$Treatment), paste(w$Block, w$Treatment)), ]
days <- sort(unique(blight_progress$Day))
chk("audpc_by_unit (16 plots) vs agricolae::audpc on the wide table", au$audpc, agricolae::audpc(w[, -(1:2)], days), 1e-12)
chk("audpc_by_unit AUDPS vs agricolae::audps (data frame)", au$audps, unlist(agricolae::audps(w[, -(1:2)], days)), 1e-12)
shuf <- blight_progress[sample(nrow(blight_progress)), ]
as <- audpc_by_unit(shuf, "Severity_pct", "Day", c("Block", "Treatment"))
chk("audpc_by_unit does not depend on the row order", as$audps[order(paste(as$Block, as$Treatment))],
    au$audps[order(paste(au$Block, au$Treatment))], 1e-12)
fa <- design_anova(au, "audps", "Treatment", block = "Block")
chk("ANOVA of AUDPS: F vs aov", fa$table$F[2], summary(aov(audps ~ factor(Block) + Treatment, au))[[1]][2, "F value"], 1e-8)

# ======================= 6. Smith's index =======================
cat("\n--- Smith's index of soil heterogeneity ---\n")
data(rice)
sm <- smith_index(rice)
ag <- index.smith(rice, PLOT = FALSE)$uniformity
N <- prod(dim(rice))
key_m <- paste(sm$table$width, sm$table$length); key_a <- paste(ag[, "Width"], ag[, "Length"])
mm <- sm$table[match(key_a, key_m), ]
chk(sprintf("rice: all %d agricolae shapes found", nrow(ag)), sum(!is.na(mm$size)), nrow(ag), 1e-12)
chk("rice: V_x = agricolae Vx (N-1)/(x (P-1))  [divisor convention]", mm$variance, ag[, "Vx"] * (N - 1) / (mm$size * (mm$plots - 1)), 1e-10)
note("agricolae divides the between-plot SS (per-unit scale) by N - 1 = %d for every shape; Smith (1938) and Gomez & Gomez use var(T/x) with P - 1 df", N - 1)
# brute force V_x of a few shapes with explicit loops over the plots
bf_v <- function(A, l, w) {
  tot <- c(); for (i in seq(1, nrow(A), by = l)) for (j in seq(1, ncol(A), by = w)) tot <- c(tot, sum(A[i:(i + l - 1), j:(j + w - 1)]))
  var(tot / (l * w))
}
R <- as.matrix(rice)
for (s in list(c(1, 1), c(2, 3), c(4, 6), c(9, 2), c(36, 1))) {
  chk(sprintf("rice: V_x of %d x %d plots vs explicit loops", s[1], s[2]),
      sm$table$variance[sm$table$length == s[1] & sm$table$width == s[2]], bf_v(R, s[1], s[2]), 1e-10)
}
u <- sm$table[sm$table$used, ]
lmw <- lm(log(variance) ~ log(size), u, weights = df)
chk("rice: b and V1 vs lm(log V ~ log x, weights = df)", c(sm$b, sm$V1), c(-coef(lmw)[2], exp(coef(lmw)[1])), 1e-10)
chk("rice: R-squared vs summary(lm)", sm$r_squared, summary(lmw)$r.squared, 1e-10)
smn <- smith_index(rice, weights = "none")
lmu <- lm(log(variance) ~ log(size), u)
chk("rice: unweighted b vs lm", smn$b, -coef(lmu)[2], 1e-10)
sms <- smith_index(rice, strips = TRUE)
chk("rice: strips = TRUE uses every shape", sum(sms$table$used), nrow(sms$table), 1e-12)
note("rice: b = %.4f with the whole-field strips (36 x 1 strips: V = %.0f) vs %.4f without", sms$b, sms$table$variance[sms$table$length == 36 & sms$table$width == 1], sm$b)
# agricolae's own b: regression through the origin of log(V1/Vx) on log x, on its Vx
za <- log(ag[-1, "Size"]); ya <- log(ag[1, "Vx"] / ag[-1, "Vx"])
b_ag <- unname(coef(lm(ya ~ 0 + za)))
note("rice: b = %.4f (weighted, strips left out) vs %.4f (agricolae's internal fit through V1 on its Vx; it does not return b)", sm$b, b_ag)
chk("rice: CV of the basic units vs agricolae (x = 1)", round(sm$cv1, 1), ag[1, "CV"], 1e-12)
# simulated independent units: b near 1
set.seed(3)
ind <- smith_index(matrix(rnorm(48 * 48, 10, 1), 48))
chk("independent units: b = 1 within 0.05", abs(ind$b - 1) < 0.05, TRUE)
# exact Smith law: V_x of a field with constant columns... (perfectly correlated units) gives b = 0 for width shapes
smc <- smith_index(uniformity_trial, cost_ratio = 3)
chk("Smith's optimum = b K1 / ((1 - b) K2)", smc$optimum, smc$b * 3 / (1 - smc$b), 1e-12)
hs <- hatheway_size(smc, reps = 4, difference = 10, treatments = 6)
t1 <- qt(0.975, 15); t2 <- qt(0.80, 15)
chk("Hatheway size vs formula", hs$plot_size, (2 * (t1 + t2)^2 * smc$cv1^2 / (4 * 10^2))^(1 / smc$b), 1e-12)

# ======================= 7. print and plot methods =======================
cat("\n--- print and plot ---\n")
out <- capture.output(print(wd), print(du), print(vw), print(mt), print(sm))
chk("print methods write their extra lines", sum(grepl("Bayesian k-ratio|Design: 7 treatments|Variance of the normal|Grand median|Fit of V_x", out)), 5, 1e-12)
pdf(NULL)
r1 <- plot(wd); r2 <- plot(du); r3 <- plot(sm); invisible(dev.off())
chk("plots return their data", c(nrow(r1), nrow(r2), nrow(r3)), c(4, 7, nrow(sm$table)), 1e-12)
chk("plot.ad_smith restores par()", par("mar"), c(5.1, 4.1, 4.1, 2.1), 1e-12)

cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
