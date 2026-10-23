# Stage 7b validation: plant-breeding analyses (AMMI, stability statistics,
# line x tester, Griffing's diallels, North Carolina designs) against
# agricolae (AMMI, index.AMMI, stability.par, stability.nonpar, lineXtester,
# carolina), base R fits (lm / anova, svd), lme4 (REML components of balanced
# data) and the closed formulas of the original papers. agricolae's data sets
# are used here only; they are not shipped with AgriDesign.
# Where the reference uses another convention or formula, the check
# confirms which formula it applies and a NOTE explains the difference.
suppressMessages({ library(agricolae); library(lme4) })
suppressMessages(devtools::load_all(quiet = TRUE))   # nc_design() does not clash with agricolae::carolina()
ok <- 0; bad <- 0
chk <- function(what, got, ref, tol = 1e-6) {
  got <- as.numeric(got); ref <- as.numeric(ref)
  d <- if (length(got) == length(ref)) max(abs(got - ref) / pmax(1, abs(ref))) else Inf
  pass <- length(got) == length(ref) && is.finite(d) && d < tol
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-66s got %-22s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(signif(utils::head(got, 4), 6), collapse = " "), paste(signif(utils::head(ref, 4), 6), collapse = " ")))
}
same <- function(what, got, ref) {
  pass <- identical(as.character(got), as.character(ref))
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-66s got %-22s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(utils::head(got, 6), collapse = " "), paste(utils::head(ref, 6), collapse = " ")))
}
note <- function(...) cat("NOTE", ..., "\n")
row_of <- function(tab, s) tab[trimws(tab$source) == s, ]
align <- function(a, b) a * sign(sum(a * b))   # the sign of an SVD axis is arbitrary

# ======================= AMMI =======================
cat("\n===== AMMI: agricolae plrv (28 genotypes x 6 localities x 3 replicates) =====\n")
data(plrv)
ag <- with(plrv, AMMI(Locality, Genotype, Rep, Yield, console = FALSE))
fit <- ammi(plrv, "Yield", "Genotype", "Locality", block = "Rep")
A <- as.data.frame(ag$ANOVA)
for (s in c("ENV", "REP(ENV)", "GEN", "GxE")) {
  rs <- if (s == "GxE") "ENV:GEN" else s
  chk(paste("AMMI ANOVA", s, ": df, SS, MS, F, p vs agricolae"),
      unlist(row_of(fit$anova, s)[c("df", "ss", "ms", "F", "p")]),
      unlist(A[rs, c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]))
}
chk("AMMI ANOVA residual: df, SS, MS vs agricolae", unlist(row_of(fit$anova, "Residuals")[c("df", "ss", "ms")]),
    unlist(A["Residuals", c("Df", "Sum Sq", "Mean Sq")]))
pl <- plrv; pl$Rep <- factor(pl$Rep)
lmA <- anova(lm(terms(Yield ~ Locality + Locality:Rep + Genotype + Genotype:Locality, keep.order = TRUE), pl))
chk("AMMI ANOVA SS vs anova(lm(env/rep + gen + gen:env))", fit$anova$ss[1:4], lmA[["Sum Sq"]][1:4])
chk("AMMI residual MS vs lm", fit$mse, lmA["Residuals", "Mean Sq"])
an <- ag$analysis
chk("AMMI axes: SS vs agricolae", fit$axes$ss, an$Sum.Sq)
chk("AMMI axes: Gollob df vs agricolae", fit$axes$df, an$Df)
chk("AMMI axes: MS vs agricolae", fit$axes$ms, an$Mean.Sq)
chk("AMMI axes: F (2 decimals) vs agricolae", round(fit$axes$F, 2), an$F.value)
chk("agricolae AMMI axes: p = pf(F rounded to 2 decimals), 4 decimals", round(pf(round(fit$axes$F, 2), fit$axes$df, fit$df_error, lower.tail = FALSE), 4), an$Pr.F)
chk("AMMI axes: p vs pf(MS / MSE) of base R", fit$axes$p, pf(an$Mean.Sq / A["Residuals", "Mean Sq"], an$Df, A["Residuals", "Df"], lower.tail = FALSE))
note("agricolae computes the p-values of the IPCA axes from the F values rounded to two decimals (IPCA5: 0.7059 instead of",
     signif(fit$axes$p[5], 4), "); AgriDesign uses the exact F.")
chk("AMMI axes: % of GxE (1 decimal) vs agricolae", round(fit$axes$pct, 1), an$percent)
bp <- ag$biplot
gb <- bp[bp$type == "GEN", ]; eb <- bp[bp$type == "ENV", ]
gb <- gb[match(fit$genotypes$genotype, trimws(rownames(gb))), ]
eb <- eb[match(fit$environments$environment, trimws(rownames(eb))), ]
chk("AMMI genotype means vs agricolae", fit$genotypes$mean, gb$Yield)
chk("AMMI environment means vs agricolae", fit$environments$mean, eb$Yield)
for (k in 1:5) {
  pc <- paste0("PC", k); ip <- paste0("IPCA", k)
  chk(paste("AMMI genotype scores", ip, "vs agricolae (sign aligned)"), align(fit$genotypes[[ip]], gb[[pc]]), gb[[pc]])
  chk(paste("AMMI environment scores", ip, "vs agricolae (sign aligned)"), align(fit$environments[[ip]], eb[[pc]]), eb[[pc]])
}
ia <- index.AMMI(ag)
ia <- ia[match(fit$genotypes$genotype, trimws(rownames(ia))), ]
chk("AMMI ASV vs agricolae index.AMMI", fit$genotypes$asv, ia$ASV)
chk("AMMI rank of ASV vs agricolae", fit$genotypes$rank_asv, ia$rASV)
chk("AMMI YSI = rank ASV + rank mean vs agricolae", fit$genotypes$ysi, ia$YSI)
# independent: svd by hand and the reconstruction of the interaction
Y <- with(plrv, tapply(Yield, list(Genotype, Locality), mean))
Z <- Y - outer(rowMeans(Y), colMeans(Y), "+") + mean(Y)
s <- svd(Z)
chk("AMMI singular values vs base svd of the interaction", fit$axes$singular_value, s$d[1:5])
G <- as.matrix(fit$genotypes[paste0("IPCA", 1:5)]); E <- as.matrix(fit$environments[paste0("IPCA", 1:5)])
chk("AMMI scores reproduce the interaction (G E')", G %*% t(E), Z[fit$genotypes$genotype, fit$environments$environment])
chk("AMMI sum of the axes = GxE sum of squares", sum(fit$axes$ss), row_of(fit$anova, "GxE")$ss)
chk("AMMI Gollob df add up to (g-1)(e-1)", sum(fit$axes$df), 27 * 5)
ASV <- sqrt((s$d[1]^2 / s$d[2]^2 * s$u[, 1] * sqrt(s$d[1]))^2 + (s$u[, 2] * sqrt(s$d[2]))^2)
chk("AMMI ASV vs the formula of Purchase et al. (2000)", fit$genotypes$asv, ASV[match(fit$genotypes$genotype, rownames(Z))])
f3 <- ammi(plrv, "Yield", "Genotype", "Locality", block = "Rep", n_axes = 2)
chk("AMMI n_axes = 2: residual GxE = IPCA3 + IPCA4 + IPCA5", row_of(f3$anova, "Residual GxE")[c("df", "ss")],
    c(sum(an$Df[3:5]), sum(an$Sum.Sq[3:5])))
# means only, with the error of the trial (agricolae's sinRepAmmi example)
data(sinRepAmmi)
if (!anyNA(sinRepAmmi$YLD) && all(table(sinRepAmmi$GEN, sinRepAmmi$ENV) == 1)) {
  agm <- with(sinRepAmmi, AMMI(ENV, GEN, 3, YLD, 93.24224, console = FALSE))
  am <- as.data.frame(agm$ANOVA)
  fm <- ammi(sinRepAmmi, "YLD", "GEN", "ENV", mse = 93.24224, df_error = am["Residuals", "Df"], reps = 3)
  chk("AMMI means only: SS GEN, GxE vs agricolae (agricolae leaves ENV empty)", fm$anova$ss[2:3], am[c("GEN", "ENV:GEN"), "Sum Sq"])
  chk("AMMI means only: F of GEN and GxE vs agricolae", fm$anova$F[2:3], am[c("GEN", "ENV:GEN"), "F value"])
  chk("AMMI means only: axes SS vs agricolae", fm$axes$ss, agm$analysis$Sum.Sq)
  chk("AMMI means only: axes F vs agricolae (2 decimals)", round(fm$axes$F, 2), agm$analysis$F.value)
} else note("sinRepAmmi has missing cells (agricolae imputes them); AMMI with means is checked on maize_met below")
mm <- aggregate(Yield_t_ha ~ Hybrid + Environment, maize_met, mean)
fb <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
fmn <- ammi(mm, "Yield_t_ha", "Hybrid", "Environment", mse = fb$mse, df_error = fb$df_error, reps = 3)
chk("AMMI means only (maize_met): ANOVA and axes equal the plot analysis",
    c(fmn$anova$ss[1:3], fmn$axes$ss, fmn$axes$p), c(fb$anova$ss[c(1, 3, 4)], fb$axes$ss, fb$axes$p))

# ======================= stability =======================
cat("\n===== Stability: agricolae dataStb (17 genotypes x 12 environments, r = 4, MSE = 1.8) =====\n")
dataStb <- read.csv(system.file("external/dataStb.csv", package = "agricolae"))
longStb <- data.frame(G = rep(rownames(dataStb), ncol(dataStb)), E = rep(names(dataStb), each = nrow(dataStb)),
                      y = unlist(dataStb), stringsAsFactors = FALSE)
gS <- nrow(dataStb); eS <- ncol(dataStb)
dfeS <- eS * 3 * (gS - 1)   # agricolae: e (r - 1)(g - 1)
sp <- stability.par(dataStb, rep = 4, MSerror = 1.8, alpha = 0.1)
st <- stability(longStb, "y", "G", "E", mse = 1.8, df_error = dfeS, reps = 4, alpha = 0.1)
an <- sp$analysis
num <- function(x) suppressWarnings(as.numeric(x))
map <- c(Genotypes = "Genotypes", Environments = "Environments", Interaction = "GxE",
         Heterogeneity = "Heterogeneity", Residual = "Residual")
for (k in names(map)) {
  r <- row_of(st$anova, map[[k]])
  chk(paste("stability ANOVA", k, ": df, SS (4 dec.), F (2 dec.) vs agricolae"),
      c(r$df, round(r$ss, 4), round(r$F, 2)), c(an[k, "Df"], num(an[k, "Sum Sq"]), num(an[k, "F value"])))
}
chk("stability ANOVA total: df, SS vs agricolae", unlist(row_of(st$anova, "Total")[c("df", "ss")]),
    c(an["Total", "Df"], num(an["Total", "Sum Sq"])))
ss <- sp$statistics
sh <- st$shukla
chk("Shukla sigma2 vs agricolae Sigma-square", sh$sigma2, ss[["Sigma-square"]])
chk("Shukla s2 (heterogeneity removed) vs agricolae s-square", sh$s2, ss[["s-square"]])
chk("Wricke ecovalence (plot scale) vs agricolae Ecovalence", st$wricke$W, ss$Ecovalence)
chk("genotype means vs agricolae", sh$mean, ss$Mean)
lab <- function(p) ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))
same("Shukla sigma2: significance (F, e-1 and error df) vs agricolae", lab(sh$p), trimws(ss[[3]]))
same("Shukla s2: significance (F, e-2 and error df) vs agricolae", lab(sh$p_s2), trimws(ss[[5]]))
ka <- sp$stability; kt <- st$kang$table
chk("Kang: rank of the mean vs agricolae", kt$rank, ka$Rank)
chk("Kang: LSD adjustment vs agricolae", kt$adjustment, ka$Adj.rank)
chk("Kang: adjusted rank vs agricolae", kt$adjusted_rank, ka$Adjusted)
chk("Kang: stability rating vs agricolae", kt$rating, ka$Stab.rating)
chk("Kang: YS vs agricolae", kt$ys, ka$YSi)
same("Kang: selected genotypes vs agricolae", ifelse(kt$selected, "+", ""), trimws(ka[["..."]]))
out <- capture.output(stability.par(dataStb, rep = 4, MSerror = 1.8, alpha = 0.1, console = TRUE))
lsd_ag <- num(sub(".*:", "", grep("LSD", out, value = TRUE)[1]))
chk("Kang: LSD vs agricolae (qt(1 - alpha/2) with alpha = 0.1)", round(st$kang$lsd, 7), lsd_ag)
note("agricolae labels this LSD '(0.05)' whatever alpha is; with alpha = 0.1 it is a two-sided 10 % LSD.",
     "AgriDesign states the alpha it used.")
# Shukla by an independent algebraic form: sigma2_i = g/((g-2)(e-1)) W_i - SS_GE/((g-1)(g-2)(e-1)), plot scale
Ym <- as.matrix(dataStb)
Zs <- Ym - outer(rowMeans(Ym), colMeans(Ym), "+") + mean(Ym)
Wi <- 4 * rowSums(Zs^2)
chk("Shukla sigma2 vs g/((g-2)(e-1)) W_i - SS_GE/((g-1)(g-2)(e-1))", sh$sigma2,
    gS / ((gS - 2) * (eS - 1)) * Wi - sum(Wi) / ((gS - 1) * (gS - 2) * (eS - 1)))
# Eberhart and Russell against lm() of every genotype on the environmental index
Ii <- colMeans(Ym) - mean(Ym)
fits <- lapply(seq_len(gS), function(i) summary(lm(Ym[i, ] ~ Ii)))
er <- st$eberhart$table
chk("Eberhart-Russell b vs lm slope", er$b, sapply(fits, function(f) f$coefficients[2, 1]))
chk("Eberhart-Russell SE(b) vs lm", er$se_b, sapply(fits, function(f) f$coefficients[2, 2]))
chk("Eberhart-Russell R2 vs lm", er$r2, sapply(fits, function(f) f$r.squared))
chk("Eberhart-Russell s2d = sigma(lm)^2 - MSE/r", er$s2d, sapply(fits, function(f) f$sigma^2) - 1.8 / 4)
chk("Eberhart-Russell t of b = 1 vs lm offset", er$t_b1,
    sapply(seq_len(gS), function(i) summary(lm(Ym[i, ] - Ii ~ Ii))$coefficients[2, 3]))
chk("Eberhart-Russell pooled deviation SS = heterogeneity residual", row_of(st$eberhart$anova, "Pooled deviation")$ss,
    row_of(st$anova, "Residual")$ss)
chk("Eberhart-Russell GxE (linear) SS vs agricolae heterogeneity", round(row_of(st$eberhart$anova, "G x E (linear)")$ss, 4),
    num(an["Heterogeneity", "Sum Sq"]))
chk("Eberhart-Russell pooled deviation SS = 4 x sum of lm residual SS", row_of(st$eberhart$anova, "Pooled deviation")$ss,
    4 * sum(sapply(fits, function(f) sum(f$residuals^2))))
xm <- colMeans(Ym)
fw <- lapply(seq_len(gS), function(i) coef(lm(Ym[i, ] ~ xm)))
chk("Finlay-Wilkinson intercepts and slopes vs lm on the environment mean", c(st$finlay$intercept, st$finlay$b),
    c(sapply(fw, `[`, 1), sapply(fw, `[`, 2)))
stl <- stability(longStb, "y", "G", "E", mse = 1.8, df_error = dfeS, reps = 4, method = "finlay", log_finlay = TRUE)
fwl <- lapply(seq_len(gS), function(i) coef(lm(log(Ym[i, ]) ~ colMeans(log(Ym)))))
chk("Finlay-Wilkinson on the log scale vs lm", stl$finlay$b, sapply(fwl, `[`, 2))
# Huehn's rank statistics
np <- stability.nonpar(cbind(name = rownames(dataStb), dataStb), "y", console = FALSE)
hr <- np$ranking; hs <- np$statistics; ht <- st$huehn$table
chk("Huehn S1 vs agricolae stability.nonpar", ht$S1, hr$s1)
chk("Huehn S2 vs agricolae", ht$S2, hr$s2)
chk("Huehn Z2 vs agricolae", ht$Z2, hr$Z2)
chk("Huehn rank of the mean vs agricolae", ht$rank_mean, hr$Rank)
chk("Huehn E(S1), E(S2), Var(S2) vs agricolae", st$huehn$expectation[c("E1", "E2", "V2")],
    unlist(hs[c("es1", "es2", "vs2")]))
K <- gS; N <- eS
chk("agricolae Var(S1) = [(K^2-4)(N+3)+30] / (45 N(N-1))  (no (K^2-1)/K^2 factor)", hs$vs1,
    ((K^2 - 4) * (N + 3) + 30) / (45 * N * (N - 1)))
chk("AgriDesign Var(S1) = (K^2-1)[(K^2-4)(N+3)+30] / (45 K^2 N(N-1))  (Nassar and Huehn 1987)", st$huehn$expectation[["V1"]],
    (K^2 - 1) * ((K^2 - 4) * (N + 3) + 30) / (45 * K^2 * N * (N - 1)))
chk("Huehn Z1 vs agricolae Z1 rescaled by the ratio of the two Var(S1)", ht$Z1,
    hr$Z1 * hs$vs1 / st$huehn$expectation[["V1"]])
# exact null distribution: the rank of a genotype is uniform on 1..K and independent across the N environments
exact <- function(K, N) {
  R <- as.matrix(expand.grid(rep(list(seq_len(K)), N)))
  s1 <- apply(R, 1, function(v) sum(abs(outer(v, v, "-"))) / (N * (N - 1)))
  s2 <- apply(R, 1, var)
  c(E1 = mean(s1), V1 = mean((s1 - mean(s1))^2), E2 = mean(s2), V2 = mean((s2 - mean(s2))^2))
}
for (kn in list(c(3, 3), c(4, 3), c(5, 4))) {
  K <- kn[1]; N <- kn[2]
  ex <- exact(K, N)
  dK <- data.frame(G = rep(seq_len(K), N), E = rep(seq_len(N), each = K), y = c(outer(seq_len(K), seq_len(N))) + 0.5)
  hK <- stability(dK, "y", "G", "E", method = "huehn")$huehn$expectation
  chk(sprintf("Huehn E(S1), Var(S1), E(S2), Var(S2) vs exact enumeration, K = %d, N = %d", K, N), hK, ex)
  cat(sprintf("      agricolae's Var(S1) would be %.6f (exact %.6f)\n", ((K^2 - 4) * (N + 3) + 30) / (45 * N * (N - 1)), ex[["V1"]]))
}
note("Huehn's Var(S1): the exact null distribution confirms Nassar and Huehn's (K^2-1)[(K^2-4)(N+3)+30]/(45K^2N(N-1));",
     "agricolae uses Var(S1) without (K^2-1)/K^2, so the two Var(S1) differ by K^2/(K^2-1) and Z1 slightly.")
st05 <- stability(longStb, "y", "G", "E", mse = 1.8, df_error = dfeS, reps = 4, method = "huehn", alpha = 0.05)
chk("Huehn critical chi-square (each Z with Bonferroni, sum) vs agricolae (alpha 0.05)", st05$huehn$critical,
    unlist(hs[c("chi.ind", "chi.sum")]), 1e-4)
# S1 and S2 by hand from the ranks of the corrected values
Rk <- apply(Ym - rowMeans(Ym) + mean(Ym), 2, rank)
S1h <- apply(Rk, 1, function(v) { s <- 0; for (j in 1:(eS - 1)) for (k in (j + 1):eS) s <- s + abs(v[j] - v[k]); 2 * s / (eS * (eS - 1)) })
chk("Huehn S1 vs double loop over pairs of environments", ht$S1, S1h)
chk("Huehn S2 vs var() of the ranks", ht$S2, apply(Rk, 1, var))
# second example (agricolae help page, 17 x 7, MSE 1.8, r = 4)
v <- list(c(10.2,8.8,8.8,9.3,9.6,7.2,8.4,9.6,7.9,10,9.3,8.0,10.1,9.4,10.8,6.3,7.4),
          c(7,7.8,7,6.9,7,8.3,7.4,6.5,6.8,7.9,7.3,6.8,8.1,7.1,7.1,6.4,4.1),
          c(5.3,4.4,5.3,4.4,5.5,4.6,6.2,6.0,6.5,5.3,5.7,4.4,4.2,5.6,5.8,3.9,3.8),
          c(7.8,5.9,7.3,5.9,7.8,6.3,7.9,7.5,7.6,5.4,5.6,7.8,6.5,8.1,7.5,5.0,5.4),
          c(9,9.2,8.8,10.6,8.3,9.3,9.6,8.8,7.9,9.1,7.7,9.5,9.4,9.4,10.3,8.8,8.7),
          c(6.9,7.7,7.9,7.9,8.0,8.1,8.2,8.6,8.8,8.9,9.2,9.4,9.5,9.6,9.9,10.2,10.4),
          c(4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.0,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7))
study <- as.data.frame(v, col.names = paste0("v", 1:7)); rownames(study) <- LETTERS[1:17]
ls2 <- data.frame(G = rep(LETTERS[1:17], 7), E = rep(names(study), each = 17), y = unlist(study))
sp2 <- stability.par(study, rep = 4, MSerror = 1.8, alpha = 0.1)
st2 <- stability(ls2, "y", "G", "E", mse = 1.8, df_error = 7 * 3 * 16, reps = 4, alpha = 0.1)
chk("study: Shukla sigma2 and s2 vs agricolae", c(st2$shukla$sigma2, st2$shukla$s2),
    c(sp2$statistics[["Sigma-square"]], sp2$statistics[["s-square"]]))
chk("study: Kang YS vs agricolae (ties in the mean)", st2$kang$table$ys, sp2$stability$YSi)
np2 <- stability.nonpar(cbind(name = LETTERS[1:17], study), "y", console = FALSE)
# exact ranks: 70 (y_ij - mean_i) are integers, so ties are exact
Ys <- as.matrix(study); Wint <- round(70 * Ys - 10 * rowSums(Ys))
Rint <- apply(Wint, 2, rank)
chk("study: Huehn S1, S2 vs ranks computed in exact integer arithmetic", c(st2$huehn$table$S1, st2$huehn$table$S2),
    c(apply(Rint, 1, function(v) sum(abs(outer(v, v, "-"))) / 42), apply(Rint, 1, var)))
Rflt <- apply(Ys - rowMeans(Ys) + mean(Ys), 2, rank)
chk("agricolae S1 = S1 from floating-point ranks (exact ties broken by rounding error)", np2$ranking$s1,
    apply(Rflt, 1, function(v) sum(abs(outer(v, v, "-"))) / 42))
note("Huehn ties: in the study data B and D are tied in environment v7 (4.4 - 6.8857 = 4.6 - 7.0857), but the corrected",
     "values differ in the 15th digit, so agricolae gives them different ranks (S1 of B: 5.0476 instead of 5.0952).",
     "AgriDesign ranks values equal to 12 significant digits as ties, which reproduces exact arithmetic.")

# ======================= line x tester =======================
cat("\n===== Line x tester: agricolae heterosis, place 1 (8 lines x 3 testers + 11 parents, r = 3) =====\n")
data(heterosis)
h1 <- subset(heterosis, Place == 1 & Factor != "Control")
out <- capture.output(agl <- with(h1, lineXtester(Replication, Female, Male, v2)))
grab <- function(pattern, from = 1) {
  i <- grep(pattern, out)
  i <- i[i >= from][1]
  x <- regmatches(out[i], gregexpr("-?[0-9]+\\.?[0-9]*(e-?[0-9]+)?", out[i]))[[1]]
  as.numeric(x)
}
after <- function(label) num(sub(".*:", "", out[grep(label, out, fixed = TRUE)[1]]))
lt <- line_tester(h1, "v2", "Female", "Male", block = "Replication")
i2 <- grep("ANOVA for line X tester analysis", out)[1]
chk("LxT ANOVA Lines: df, SS, MS, F (vs LxT) vs agricolae", unlist(row_of(lt$anova_lxt, "Lines")[c("df", "ss", "ms")]),
    grab("^Lines ", i2)[1:3])
chk("LxT Lines F against L x T (3 decimals) vs agricolae", round(row_of(lt$anova_lxt, "Lines")$F, 3), grab("^Lines ", i2)[4])
chk("LxT ANOVA Testers: df, SS, MS vs agricolae", unlist(row_of(lt$anova_lxt, "Testers")[c("df", "ss", "ms")]),
    grab("^Testers ", i2)[1:3])
chk("LxT Testers F against L x T vs agricolae", round(row_of(lt$anova_lxt, "Testers")$F, 3), grab("^Testers ", i2)[4])
chk("LxT ANOVA L x T: df, SS, MS, F vs agricolae", unlist(row_of(lt$anova_lxt, "Lines x Testers")[c("df", "ss", "ms")]),
    grab("^Lines X Testers", i2)[1:3])
for (s in c("Replications", "Treatments", "Parents", "Parents vs. Crosses", "Crosses", "Error", "Total")) {
  mine <- c(Replications = "Blocks", Treatments = "Entries", Parents = "Parents", "Parents vs. Crosses" = "Parents vs crosses",
            Crosses = "Crosses", Error = "Error", Total = "Total")[[s]]
  chk(paste("LxT ANOVA", s, ": df, SS vs agricolae"), unlist(row_of(lt$anova, mine)[c("df", "ss")]),
      grab(paste0("^", gsub(".", "\\.", s, fixed = TRUE), " "))[1:2])
}
h1f <- h1; h1f$Replication <- factor(h1f$Replication)
h1f$entry <- factor(paste(h1f$Female, h1f$Male))
mE <- lm(v2 ~ Replication + entry, h1f)
chk("LxT error MS vs lm(y ~ block + entry)", lt$mse, deviance(mE) / df.residual(mE))
cr <- subset(h1f, !is.na(Female) & !is.na(Male))
lmc <- anova(lm(v2 ~ Replication + Female * Male, cr))
chk("LxT lines, testers, LxT SS vs anova(lm) of the crosses", c(row_of(lt$anova, "Lines")$ss, row_of(lt$anova, "Testers")$ss,
                                                                 row_of(lt$anova, "Lines x Testers")$ss), lmc[2:4, "Sum Sq"])
il <- grep("Lines Effects:", out); it <- grep("Testers Effects:", out)
nmL <- strsplit(trimws(out[il + 1]), "\\s+")[[1]]; vL <- as.numeric(strsplit(trimws(out[il + 2]), "\\s+")[[1]])
nmT <- strsplit(trimws(out[it + 1]), "\\s+")[[1]]; vT <- as.numeric(strsplit(trimws(out[it + 2]), "\\s+")[[1]])
chk("LxT GCA of lines (3 decimals) vs agricolae", round(lt$gca_lines$gca[match(nmL, lt$gca_lines$line)], 3), vL)
chk("LxT GCA of testers (3 decimals) vs agricolae", round(lt$gca_testers$gca[match(nmT, lt$gca_testers$tester)], 3), vT)
is_ <- grep("SCA Effects:", out)
hdr <- strsplit(trimws(out[is_ + 3]), "\\s+")[[1]][-1]
scag <- do.call(rbind, lapply(out[is_ + 3 + seq_along(nmL)], function(z) {
  p <- strsplit(trimws(z), "\\s+")[[1]]; data.frame(line = p[1], tester = hdr, sca = as.numeric(p[-1]))
}))
mine <- lt$sca[match(paste(scag$line, scag$tester), paste(lt$sca$line, lt$sca$tester)), ]
chk("LxT SCA effects (3 decimals) vs agricolae", round(mine$sca, 3), scag$sca)
chk("LxT SE of GCA line, GCA tester, SCA vs agricolae", lt$se[1:3],
    c(after("S.E. (gca for line)"), after("S.E. (gca for tester)"), after("S.E. (sca effect)")))
chk("LxT SE of differences (lines, testers, SCA) vs agricolae", lt$se[4:6],
    c(after("S.E. (gi - gj)line"), after("S.E. (gi - gj)tester"), after("S.E. (sij - skl)tester")))
chk("LxT Cov HS (line), Cov HS (tester) vs agricolae", lt$components$estimate[1:2],
    c(after("Cov H.S. (line)"), after("Cov H.S. (tester)")), 1e-5)
chk("LxT contributions of lines, testers, LxT vs agricolae", lt$contribution,
    c(after("Contributions of lines"), after("Contributions of testers"), after("Contributions of lxt")), 1e-5)
# the average covariance of half sibs and what follows from it
l <- 8; t <- 3; r <- 3
Ml <- row_of(lt$anova_lxt, "Lines")$ms; Mt <- row_of(lt$anova_lxt, "Testers")$ms
Mlt <- row_of(lt$anova_lxt, "Lines x Testers")$ms; Me <- lt$mse
wav <- ((l - 1) * Ml + (t - 1) * Mt) / (l + t - 2)
chk("agricolae Cov HS (average) = [wavg - MSlt] / (r(2lt - l - t))  (Singh-Chaudhary as printed)",
    (wav - Mlt) / (r * (2 * l * t - l - t)), after("Cov H.S. (average)"), 1e-5)
chk("AgriDesign Cov HS (average) = [wavg - MSlt] (l+t-2) / (r(2lt - l - t))  (moment estimator)",
    lt$components$estimate[3], (wav - Mlt) * (l + t - 2) / (r * (2 * l * t - l - t)))
chk("AgriDesign Cov HS (average) = d.f.-weighted mean of the two moment estimators",
    lt$components$estimate[3], ((l - 1) * t * lt$components$estimate[1] + (t - 1) * l * lt$components$estimate[2]) /
      ((l - 1) * t + (t - 1) * l))
note("agricolae's Cov HS (average) uses the weighted mean square without the factor (l + t - 2): it gives",
     signif(after("Cov H.S. (average)"), 5), "although Cov HS of lines and testers are", signif(lt$components$estimate[1], 5),
     "and", signif(lt$components$estimate[2], 5), "; E[((l-1)MSl + (t-1)MSt)/(l+t-2)] = Ve + r Vlt + r (2lt-l-t)/(l+t-2) Cov HS,",
     "so AgriDesign gives", signif(lt$components$estimate[3], 5), "(it lies between the two, as the expectations require).")
chk("agricolae Cov FS = sum of (MS - MSe) / 3r, without the Cov HS terms",
    ((Ml - Me) + (Mt - Me) + (Mlt - Me)) / (3 * r), after("Cov F.S. (average)"), 1e-5)
chk("AgriDesign Cov FS = [sum (MS - MSe) + 6r CovHS - r(l+t) CovHS] / 3r (Singh-Chaudhary)", lt$components$estimate[4],
    ((Ml - Me) + (Mt - Me) + (Mlt - Me) + (6 * r - r * (l + t)) * lt$components$estimate[3]) / (3 * r))
note("agricolae's Cov FS does not include the term (6r - r(l+t)) Cov HS of Singh and Chaudhary's formula; with sigma2_l = sigma2_t = Cov HS",
     "the expectations give Cov FS = 2 Cov HS + sigma2_lt, which AgriDesign reproduces.")
chk("agricolae sigma2_A (F = 0, 1) = 4 Cov HS / (1 + F) with its own Cov HS", c(4, 2) * after("Cov H.S. (average)"),
    c(after("F = 0, Adittive genetic variance"), after("F = 1, Adittive genetic variance")), 1e-5)
s2sca <- (Mlt - Me) / r
chk("agricolae sigma2_D (F = 0, 1) = 2 sigma2_SCA / (1 + F)", c(2, 1) * s2sca,
    c(after("F = 0, Variance due to Dominance"), after("F = 1, Variance due to Dominance")), 1e-5)
chk("AgriDesign sigma2_D (F = 0, 1) = 4 sigma2_SCA / (1 + F)^2", lt$genetic$sigma2_D, c(4, 1) * s2sca)
note("Cov FS - 2 Cov HS = sigma2_SCA = ((1 + F)/2)^2 sigma2_D, so sigma2_D = 4 sigma2_SCA/(1 + F)^2: 4 sigma2_SCA",
     "for non-inbred parents. agricolae uses 2/(1 + F) without squaring, half the right value at F = 0 (equal at F = 1).")
# REML of the crosses alone (balanced and positive components: REML = moment estimators)
ltc <- line_tester(cr, "v2", "Female", "Male", block = "Replication")
cr$L <- factor(cr$Female); cr$T <- factor(cr$Male)
mr <- lmer(v2 ~ Replication + (1 | L) + (1 | T) + (1 | L:T), cr, REML = TRUE)
vc <- as.data.frame(VarCorr(mr)); vv <- setNames(vc$vcov, vc$grp)
if (all(ltc$components$estimate[c(1, 2, 6)] > 0)) {
  chk("LxT Cov HS lines, testers, sigma2 SCA vs lme4 REML (crosses only)", ltc$components$estimate[c(1, 2, 6)],
      vv[c("L", "T", "L:T")], 1e-4)
  chk("LxT error MS vs lme4 residual variance (crosses only)", ltc$mse, vv[["Residual"]], 1e-4)
} else note("a component is negative: REML truncates it at zero, comparison skipped")

# ======================= Griffing's diallels =======================
cat("\n===== Diallel (Griffing 1956), maize_diallel: 6 parents, r = 3, all four methods =====\n")
dd <- maize_diallel
P <- sort(unique(dd$Female)); p <- length(P)
cellm <- with(dd, tapply(Yield_t_ha, list(factor(Female, P), factor(Male, P)), mean))  # x[i, j]: female i, male j
# Griffing's closed formulas on the p x p table of cell means x (female in rows); cells outside the method are ignored
closed <- function(x, meth) {
  if (meth == 1) {
    Xi <- rowSums(x) + colSums(x); X <- sum(x)
    list(SSg = sum(Xi^2) / (2 * p) - 2 * X^2 / p^2,
         SSs = sum(x * (x + t(x))) / 2 - sum(Xi^2) / (2 * p) + X^2 / p^2,
         SSr = sum(((x - t(x))[upper.tri(x)])^2) / 2,
         g = Xi / (2 * p) - X / p^2, s = (x + t(x)) / 2 - outer(Xi, Xi, "+") / (2 * p) + X / p^2)
  } else if (meth == 2) {
    xu <- x; xu[lower.tri(xu)] <- 0; xs <- xu + t(xu); diag(xs) <- diag(x)
    Xi <- rowSums(xs); X <- sum(xu); xii <- diag(x)
    list(SSg = (sum((Xi + xii)^2) - 4 * X^2 / p) / (p + 2),
         SSs = sum(xu^2) - sum((Xi + xii)^2) / (p + 2) + 2 * X^2 / ((p + 1) * (p + 2)),
         g = (Xi + xii - 2 * X / p) / (p + 2), s = xs - outer(Xi + xii, Xi + xii, "+") / (p + 2) + 2 * X / ((p + 1) * (p + 2)))
  } else if (meth == 3) {
    xo <- x; diag(xo) <- 0
    Xi <- rowSums(xo) + colSums(xo); X <- sum(xo)
    list(SSg = sum(Xi^2) / (2 * (p - 2)) - 2 * X^2 / (p * (p - 2)),
         SSs = sum(((xo + t(xo))[upper.tri(xo)])^2) / 2 - sum(Xi^2) / (2 * (p - 2)) + X^2 / ((p - 1) * (p - 2)),
         SSr = sum(((xo - t(xo))[upper.tri(xo)])^2) / 2,
         g = (p * Xi - 2 * X) / (2 * p * (p - 2)), s = (xo + t(xo)) / 2 - outer(Xi, Xi, "+") / (2 * (p - 2)) + X / ((p - 1) * (p - 2)))
  } else {
    xu <- x; xu[!upper.tri(xu)] <- 0; xs <- xu + t(xu)
    Xi <- rowSums(xs); X <- sum(xu)
    list(SSg = sum(Xi^2) / (p - 2) - 4 * X^2 / (p * (p - 2)),
         SSs = sum(xu^2) - sum(Xi^2) / (p - 2) + 2 * X^2 / ((p - 1) * (p - 2)),
         g = (p * Xi - 2 * X) / (p * (p - 2)), s = xs - outer(Xi, Xi, "+") / (p - 2) + 2 * X / ((p - 1) * (p - 2)))
  }
}
# variances of Griffing (1956) as tabulated by Singh and Chaudhary (1985), in units of the error of an entry mean
textbook <- list(
  list(g = (p - 1) / (2 * p^2), sij = (p^2 - 2 * p + 2) / (2 * p^2), sii = (p - 1)^2 / p^2),
  list(g = (p - 1) / (p * (p + 2)), sij = (p^2 + p + 2) / ((p + 1) * (p + 2)), sii = p * (p - 1) / ((p + 1) * (p + 2))),
  list(g = (p - 1) / (2 * p * (p - 2)), sij = (p - 3) / (2 * (p - 1))),
  list(g = (p - 1) / (p * (p - 2)), sij = (p - 3) / (p - 1)))
for (meth in 1:4) {
  dl <- diallel(dd, "Yield_t_ha", "Female", "Male", block = "Block", method = meth)
  sel <- switch(meth, dd, subset(dd, Female <= Male), subset(dd, Female != Male), subset(dd, Female < Male))
  sel$E <- factor(paste(sel$Female, sel$Male)); sel$B <- factor(sel$Block)
  av <- anova(lm(Yield_t_ha ~ B + E, sel))
  mse <- av["Residuals", "Mean Sq"]; s2 <- mse / 3
  chk(sprintf("diallel m%d: entries and error vs anova(lm(y ~ block + entry))", meth),
      c(row_of(dl$anova_entries, "Entries")$ss, dl$mse, dl$df_error), c(av["E", "Sum Sq"], mse, av["Residuals", "Df"]))
  ss <- setNames(dl$anova$ss, trimws(dl$anova$source))
  cf <- closed(cellm, meth)
  ut <- if (meth %in% 1:2) upper.tri(cf$s, diag = TRUE) else upper.tri(cf$s)
  sca <- matrix(NA, p, p); sca[cbind(match(dl$sca$parent1, P), match(dl$sca$parent2, P))] <- dl$sca$estimate
  chk(sprintf("diallel m%d: SS GCA, SCA vs Griffing's closed formulas", meth), ss[c("GCA", "SCA")], c(cf$SSg, cf$SSs))
  chk(sprintf("diallel m%d: GCA effects vs closed formula", meth), dl$gca$estimate, cf$g)
  chk(sprintf("diallel m%d: SCA effects vs closed formula", meth), sca[ut], cf$s[ut])
  if (meth %in% c(1, 3)) {
    chk(sprintf("diallel m%d: SS reciprocal vs (1/2) sum (x_ij - x_ji)^2", meth), ss[["Reciprocal"]], cf$SSr)
    chk(sprintf("diallel m%d: maternal + non-maternal = reciprocal", meth), ss[["Maternal"]] + ss[["Non-maternal"]], ss[["Reciprocal"]])
    xo <- cellm; diag(xo) <- 0
    chk(sprintf("diallel m%d: maternal effects m_i = (X_i. - X_.i) / 2p", meth), dl$maternal$estimate,
        (rowSums(xo) - colSums(xo)) / (2 * p))
    rr <- (cellm - t(cellm)) / 2
    chk(sprintf("diallel m%d: reciprocal effects r_ij = (x_ij - x_ji)/2", meth), dl$reciprocal$estimate,
        rr[cbind(match(dl$reciprocal$parent1, P), match(dl$reciprocal$parent2, P))])
    chk(sprintf("diallel m%d: SE(r_ij) = sqrt(s2 / 2), SE(r_ij - r_kl) = sqrt(s2)", meth),
        c(dl$reciprocal$se[1], dl$se_differences$se[dl$se_differences$difference == "r_ij - r_kl"]), sqrt(c(s2 / 2, s2)))
  }
  # exact variances of the closed formulas: they are linear in the independent cell means (variance s2 each)
  used <- switch(meth, matrix(TRUE, p, p), upper.tri(cellm, diag = TRUE), row(cellm) != col(cellm), upper.tri(cellm))
  J <- sapply(which(used), function(k) { z <- matrix(0, p, p); z[k] <- 1; cz <- closed(z, meth); c(cz$g, cz$s[ut]) })
  V <- s2 * J %*% t(J)
  idx <- which(ut, arr.ind = TRUE); sp <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  ord <- order(idx[, 1], idx[, 2])   # the package lists the pairs row by row
  vg <- diag(V)[1:p]; vs <- diag(V)[p + ord]
  off <- sp[, 1] != sp[, 2]
  chk(sprintf("diallel m%d: SE of GCA effects vs exact variance of the closed formula", meth), dl$gca$se, sqrt(vg))
  chk(sprintf("diallel m%d: SE of SCA effects vs exact variance of the closed formula", meth), dl$sca$se, sqrt(vs))
  tb <- textbook[[meth]]
  chk(sprintf("diallel m%d: Var(g_i), Var(s_ij) vs Singh and Chaudhary's table", meth), c(vg[1], vs[off][1]) / s2, c(tb$g, tb$sij))
  if (!is.null(tb$sii)) chk(sprintf("diallel m%d: Var(s_ii) vs Singh and Chaudhary's table", meth), vs[!off][1] / s2, tb$sii)
  vd <- function(a, b) V[a, a] + V[b, b] - 2 * V[a, b]
  pr <- t(combn(p, 2)); dg <- mean(mapply(vd, pr[, 1], pr[, 2]))
  sv <- p + ord; np_ <- nrow(sp); pp <- t(combn(np_, 2))
  share <- apply(pp, 1, function(z) length(intersect(sp[z[1], ], sp[z[2], ])) > 0)
  s1 <- !off[pp[, 1]]; s2_ <- !off[pp[, 2]]
  grp <- list(s1 & s2_, !s1 & !s2_ & share, !s1 & !s2_ & !share)
  ref <- c(dg, vapply(grp, function(g) if (any(g)) mean(mapply(vd, sv[pp[g, 1]], sv[pp[g, 2]])) else NA_real_, 0))
  ref <- ref[!is.na(ref)]
  got <- dl$se_differences$se[dl$se_differences$difference != "r_ij - r_kl"]
  chk(sprintf("diallel m%d: SE of differences vs exact variances of the closed formulas", meth), got, sqrt(ref))
  cat(sprintf("      variances of differences / s2: %s\n", paste(signif(ref / s2, 5), collapse = ", ")))
  # model II with Griffing's expected mean squares
  ms <- setNames(dl$anova$ms, trimws(dl$anova$source)); me <- s2
  if (meth == 1) {
    s2s <- (ms[["SCA"]] - me) * p^2 / (2 * (p^2 - p + 1))
    s2g <- (ms[["GCA"]] - me - 2 * (p - 1) / p * s2s) / (2 * p)
  } else {
    s2s <- (ms[["SCA"]] - me) / c(1, 2, 1)[meth - 1]
    s2g <- (ms[["GCA"]] - ms[["SCA"]]) / c(p + 2, 2 * (p - 2), p - 2)[meth - 1]
  }
  chk(sprintf("diallel m%d: Model II sigma2 GCA, SCA from the expected mean squares", meth), dl$model2$estimate[1:2], c(s2g, s2s))
  chk(sprintf("diallel m%d: GCA vs SCA (Model II F)", meth), dl$gca_model2$F, ms[["GCA"]] / ms[["SCA"]])
  # independent check with lm(): GCA after the mean, SCA after GCA, on the cell means
  cm <- dl$means; ci <- match(cm$parent1, P); cj <- match(cm$parent2, P)
  Gm <- sapply(seq_len(p), function(k) (ci == k) + (cj == k))
  pairk <- factor(paste(pmin(ci, cj), pmax(ci, cj)))
  lmg <- suppressWarnings(anova(lm(cm$mean ~ Gm + pairk)))
  chk(sprintf("diallel m%d: SS GCA and SCA vs anova(lm(mean ~ G + pair))", meth), ss[c("GCA", "SCA")], lmg[1:2, "Sum Sq"])
}

# ======================= North Carolina designs =======================
cat("\n===== North Carolina designs vs agricolae carolina (data DC) =====\n")
data(DC)
c1 <- DC$carolina1
invisible(capture.output(ag1 <- agricolae::carolina(model = 1, c1)))
nc1 <- nc_design(c1, "yield", "male", "female", block = "rep", set = "set", design = 1)
a1 <- anova(ag1$model)
chk("NC I: SS set, rep(set), male(set), female(male), plot error, within vs agricolae", nc1$anova$ss[1:6], a1[["Sum Sq"]])
chk("NC I: df vs agricolae", nc1$anova$df[1:6], a1$Df)
chk("NC I: Vm, Vf vs agricolae var.m, var.f", nc1$components$estimate, c(ag1$var.m, ag1$var.f))
chk("NC I: sigma2_A, sigma2_D vs agricolae var.A, var.D", nc1$genetic$estimate[3:4], c(ag1$var.A, ag1$var.D))
note("NC I: agricolae tests every line against the within-plot residual (",
     signif(a1["set:male", "F value"], 5), "for males); AgriDesign tests males against females within males (F =",
     signif(row_of(nc1$anova, "Males(sets)")$F, 5), ") and females against the plot error, as the expected mean squares require.")
chk("NC I: F males = MS males / MS females(males)", row_of(nc1$anova, "Males(sets)")$F,
    a1["set:male", "Mean Sq"] / a1["set:male:female", "Mean Sq"])
chk("NC I: F females = MS females / MS plot error", row_of(nc1$anova, "Females(males)")$F,
    a1["set:male:female", "Mean Sq"] / a1["set:replication:male:female", "Mean Sq"])
c2 <- subset(DC$carolina2, Loc == 1)
invisible(capture.output(ag2 <- agricolae::carolina(model = 2, c2[, c("set", "male", "female", "rep", "yield")])))
nc2 <- nc_design(c2, "yield", "male", "female", block = "rep", set = "set", design = 2)
a2 <- anova(ag2$model)
chk("NC II: SS set, rep(set), male, female, male x female, error vs agricolae", nc2$anova$ss[1:6], a2[["Sum Sq"]])
chk("NC II: df vs agricolae", nc2$anova$df[1:6], a2$Df)
mf2 <- table(unique(c2[c("set", "female")])$set)[[1]]; mm2 <- table(unique(c2[c("set", "male")])$set)[[1]]
chk("NC II: Vm = (MSm - MSmf)/(r f), Vf, Vmf vs the expected mean squares",
    nc2$components$estimate,
    c((a2["set:male", "Mean Sq"] - a2["set:male:female", "Mean Sq"]) / (3 * mf2),
      (a2["set:female", "Mean Sq"] - a2["set:male:female", "Mean Sq"]) / (3 * mm2),
      (a2["set:male:female", "Mean Sq"] - a2["Residuals", "Mean Sq"]) / 3))
chk("NC II: Vmf vs agricolae var.mf", nc2$components$estimate[3], ag2$var.mf)
chk("NC II: sigma2_D = 4 Vmf vs agricolae var.D", nc2$genetic$estimate[4], ag2$var.D)
chk("agricolae NC II var.m = (MSm - MSmf) / (r x total number of females)", ag2$var.m,
    (a2["set:male", "Mean Sq"] - a2["set:male:female", "Mean Sq"]) / (3 * length(unique(c2$female))))
chk("agricolae NC II var.f = (MSf - MSmf) / (r x total number of males)", ag2$var.f,
    (a2["set:female", "Mean Sq"] - a2["set:male:female", "Mean Sq"]) / (3 * length(unique(c2$male))))
note("NC II: E(MS males within sets) = Ve + r Vmf + r f Vm with f = females per set (", mf2, "); agricolae divides by the",
     "number of female labels in the whole data (", length(unique(c2$female)), "); the two agree with a single set and",
     "differ by the number of sets when there are several.")
c3 <- DC$carolina3
invisible(capture.output(ag3 <- agricolae::carolina(model = 3, c3)))
nc3 <- nc_design(c3, "yield", "male", "female", block = "rep", set = "set", design = 3)
a3 <- anova(ag3$model)
mine3 <- setNames(nc3$anova$ss, nc3$anova$source)
chk("NC III: SS set, rep(set), testers within sets, males(set), males x testers, error vs agricolae",
    c(mine3[["Sets"]], mine3[["Blocks(sets)"]], mine3[["Testers"]] + mine3[["Testers x Sets"]], mine3[["Males(sets)"]],
      mine3[["Males x Testers(sets)"]], mine3[["Error"]]), a3[["Sum Sq"]])
chk("NC III: Vm = (MSm - MSe)/(2r), Vmt = (MSmt - MSe)/r vs agricolae var.m, var.mi", nc3$components$estimate,
    c(ag3$var.m, ag3$var.mi))
chk("NC III: sigma2_A = 4 Vm vs agricolae var.A", nc3$genetic$estimate[1], ag3$var.A)
chk("agricolae NC III var.D = 2 Vmt", ag3$var.D, 2 * nc3$components$estimate[2])
chk("AgriDesign NC III sigma2_D = Vmt", nc3$genetic$estimate[2], nc3$components$estimate[2])
note("NC III: with an F2 male crossed to both parents, one locus gives backcross sums (a+d, d, d-a) and differences",
     "(a-d, a, a+d) with frequencies 1/4, 1/2, 1/4: Var(sum) = a^2/2 = sigma2_A and Var(difference) = d^2/2 = 2 sigma2_D.",
     "The interaction mean square has expectation Ve + r Var(difference)/2 = Ve + r sigma2_D, so sigma2_D = (MSmt - MSe)/r = Vmt.",
     "agricolae reports var.D = 2 var.mi, which corresponds to defining Vmt as the variance of one interaction effect,",
     "with coefficient 2r instead of r; under the model above the two estimates differ by a factor of 2.")
# design III on an exactly enumerated F2: two loci, the 16 gamete combinations (frequencies 1:2:1 at each locus),
# no error; the sample variances (divisor n - 1) are 16/15 of the population ones
aa <- c(0.30, 0.20); dd3 <- c(0.25, -0.10)
G3 <- as.matrix(expand.grid(l1 = c(0, 1, 1, 2), l2 = c(0, 1, 1, 2)))
fam1 <- (G3 / 2) %*% aa + (1 - G3 / 2) %*% dd3; fam2 <- (G3 / 2) %*% dd3 - (1 - G3 / 2) %*% aa
sim <- expand.grid(Tester = c("P1", "P2"), Male = seq_len(16), Rep = 1:2, stringsAsFactors = FALSE)
sim$y <- ifelse(sim$Tester == "P1", fam1[sim$Male], fam2[sim$Male]) + c(0.1, -0.1)[sim$Rep]
n3 <- nc_design(sim, "y", "Male", "Tester", block = "Rep", design = 3)
chk("NC III exact F2 (2 loci): sigma2_A x 15/16 = sum(a^2)/2", n3$genetic$estimate[1] * 15 / 16, sum(aa^2) / 2)
chk("NC III exact F2: sigma2_D x 15/16 = sum(d^2)/4 (2 Vmt would give twice this value)", n3$genetic$estimate[2] * 15 / 16,
    sum(dd3^2) / 4)
# REML of balanced data equals the moment estimators when every component is positive
for (k in 1:3) {
  dk <- nc_designs[[k]]; dk$S <- factor(dk$Set); dk$R <- factor(paste(dk$Set, dk$Rep))
  dk$M <- factor(paste(dk$Set, dk$Male))
  if (k == 1) {
    fk <- nc_design(dk, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 1)
    dk$F <- factor(paste(dk$Set, dk$Male, dk$Female))
    mr <- lmer(Yield_t_ha ~ S + (1 | R) + (1 | M) + (1 | F), dk); ref <- c("M", "F")
  } else if (k == 2) {
    fk <- nc_design(dk, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 2)
    dk$F <- factor(paste(dk$Set, dk$Female)); dk$MF <- factor(paste(dk$Set, dk$Male, dk$Female))
    mr <- lmer(Yield_t_ha ~ S + (1 | R) + (1 | M) + (1 | F) + (1 | MF), dk); ref <- c("M", "F", "MF")
  } else next
  vc <- as.data.frame(VarCorr(mr)); vv <- setNames(vc$vcov, vc$grp)
  if (all(fk$components$estimate > 0) && vv[["R"]] > 0)
    chk(sprintf("NC %s (nc_designs): components vs lme4 REML", c("I", "II")[k]), fk$components$estimate, vv[ref], 1e-4)
  else note("NC", k, ": a component is not positive, REML comparison skipped")
}
dk <- nc_designs$nc3
f3 <- nc_design(dk, "Yield_t_ha", "Male", "Tester", block = "Rep", set = "Set", design = 3)
dk$S <- factor(dk$Set); dk$R <- factor(paste(dk$Set, dk$Rep)); dk$M <- factor(paste(dk$Set, dk$Male))
# design III with testers fixed: the sums and differences of every plant (restricted model)
pm <- with(dk, tapply(Yield_t_ha, list(M, Tester), mean)); sm <- pm[, 1] + pm[, 2]; df_ <- pm[, 1] - pm[, 2]
setm <- sub(" .*", "", rownames(pm))
vS <- sum(tapply(sm, setm, function(v) sum((v - mean(v))^2))) / (length(sm) - 2)
vD <- sum(tapply(df_, setm, function(v) sum((v - mean(v))^2))) / (length(sm) - 2)
mse3 <- row_of(f3$anova, "Error")$ms
chk("NC III: sigma2_A = Var(sum within sets) - 2 MSe/r", f3$genetic$estimate[1], vS - 2 * mse3 / 2)
chk("NC III: sigma2_D = [Var(difference within sets) - 2 MSe/r] / 2", f3$genetic$estimate[2], (vD - 2 * mse3 / 2) / 2)

cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
