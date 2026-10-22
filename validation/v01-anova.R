# Stage 1 validation: design_anova() against base R aov() with Error() strata,
# car::Anova(type = 3), agricolae and the figures of the web application.
suppressMessages(devtools::load_all(quiet = TRUE))
ok <- 0; bad <- 0
chk <- function(what, got, ref, tol = 1e-6) {
  d <- max(abs(got - ref) / pmax(1, abs(ref)))
  pass <- is.finite(d) && d < tol
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %-55s got %-28s ref %s\n", if (pass) "OK" else "FAIL", what,
              paste(signif(got, 7), collapse = " "), paste(signif(ref, 7), collapse = " ")))
}
fac <- function(d, v) { for (x in v) d[[x]] <- factor(d[[x]]); d }
row_of <- function(fit, term, col) fit$table[fit$table$term == term, col]

# --- RCBD maize: app F(3,9) = 437.267, CV 2.1 %, means 4.215 5.845 6.885 7.183
f <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block", design = "rcbd")
d <- fac(maize_nitrogen, c("Block", "Nitrogen_kg_ha"))
a <- summary(aov(Yield_t_ha ~ Block + Nitrogen_kg_ha, d))[[1]]
chk("RCBD maize: SS N vs aov", row_of(f, "Nitrogen_kg_ha", "ss"), a[2, "Sum Sq"])
chk("RCBD maize: F N vs app 437.267", round(row_of(f, "Nitrogen_kg_ha", "F"), 3), 437.267)
chk("RCBD maize: CV vs app 2.1", round(f$cv, 1), 2.1)
m <- ls_means(f, "Nitrogen_kg_ha")
chk("RCBD maize: LS means vs app", round(m$mean, 3), c(4.215, 5.845, 6.885, 7.183), 1e-3)
chk("RCBD maize: SE mean = sqrt(MSE/4)", m$se, rep(sqrt(f$residual[["ms"]] / 4), 4))

# --- CRD beans
f <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
a <- summary(aov(Yield_kg_plot ~ Variety, bean_varieties))[[1]]
chk("CRD bean: F vs aov", row_of(f, "Variety", "F"), a[1, "F value"])
chk("CRD bean: residual MS vs aov", f$residual[["ms"]], a[2, "Mean Sq"])

# --- Latin square: app F(3,6) = 164.10
f <- design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column")
d <- fac(wheat_latin, c("Row", "Column"))
a <- summary(aov(Yield_t_ha ~ Row + Column + Treatment, d))[[1]]
chk("Latin: suggested design", as.numeric(f$design == "latin"), 1)
chk("Latin: F treatment vs aov", row_of(f, "Treatment", "F"), a[3, "F value"])
chk("Latin: F vs app 164.10", round(row_of(f, "Treatment", "F"), 2), 164.10)

# --- Factorial in RCBD (tomato), Type III vs car::Anova on sum-to-zero coding
f <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block")
d <- fac(tomato_factorial, c("Block", "Variety", "Fertilizer"))
lmf <- lm(Yield_kg_plant ~ Block + Variety * Fertilizer, d,
          contrasts = list(Block = contr.sum, Variety = contr.sum, Fertilizer = contr.sum))
ca <- car::Anova(lmf, type = 3)
chk("Factorial: SS III vs car::Anova", f$table$ss, ca[c("Block", "Variety", "Fertilizer", "Variety:Fertilizer"), "Sum Sq"])
a <- summary(aov(Yield_kg_plant ~ Block + Variety * Fertilizer, d))[[1]]
f1 <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block", ss_type = 1)
chk("Factorial: SS I vs aov", f1$table$ss, a[1:4, "Sum Sq"])

# --- Split-plot in RCBD: app riego F(1,2) = 10460.6
f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd")
d <- fac(split_irrigation, "Block")
a <- summary(aov(Yield_t_ha ~ Block + Irrigation * Variety + Error(Block:Irrigation), d))
chk("Split RCBD: F Irrigation vs aov Error()", row_of(f, "Irrigation", "F"),
    a[["Error: Block:Irrigation"]][[1]]["Irrigation", "F value"])
chk("Split RCBD: F Irrigation vs app 10460.6", round(row_of(f, "Irrigation", "F"), 1), 10460.6)
w <- a[["Error: Within"]][[1]]
chk("Split RCBD: F Variety, A:B vs aov", c(row_of(f, "Variety", "F"), row_of(f, "Irrigation:Variety", "F")),
    w[c("Variety         ", "Irrigation:Variety"), "F value"])
sp <- with(d, agricolae::sp.plot(Block, Irrigation, Variety, Yield_t_ha))
chk("Split RCBD: Error a, Error b vs agricolae::sp.plot",
    c(row_of(f, "Block:Irrigation", "ms"), f$residual[["ms"]]), c(sp$gl.a, sp$gl.b) * 0 + c(sp$Ea, sp$Eb))

# --- Split-plot in CRD: Error a must be replicates within main plots: a(r-1) df
f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_crd")
a <- summary(aov(Yield_t_ha ~ Irrigation * Variety + Error(Irrigation:Block), d))
chk("Split CRD: df Error a = a(r-1)", row_of(f, "Irrigation:Block", "df"),
    a[["Error: Irrigation:Block"]][[1]]["Residuals", "Df"])
chk("Split CRD: F Irrigation vs aov Error()", row_of(f, "Irrigation", "F"),
    a[["Error: Irrigation:Block"]][[1]]["Irrigation", "F value"])

# --- Strip-plot: three strata vs aov and agricolae::strip.plot
f <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "strip_rcbd")
a <- summary(aov(Yield_t_ha ~ Block + Irrigation * Variety + Error(Block:Irrigation + Block:Variety), d))
chk("Strip: F Irrigation vs aov", row_of(f, "Irrigation", "F"), a[["Error: Block:Irrigation"]][[1]]["Irrigation", "F value"])
chk("Strip: F Variety vs aov", row_of(f, "Variety", "F"), a[["Error: Block:Variety"]][[1]]["Variety", "F value"])
chk("Strip: F A:B vs aov", row_of(f, "Irrigation:Variety", "F"), a[["Error: Within"]][[1]]["Irrigation:Variety", "F value"])

# --- Nested: plants within plots (inner labels repeat within each outer level)
set.seed(1)
nd <- expand.grid(rep = 1:3, Plant = c("p1", "p2", "p3"), Orchard = c("O1", "O2", "O3", "O4"))
nd$y <- rnorm(nrow(nd), 10 + as.integer(nd$Orchard)) + rnorm(12)[as.integer(interaction(nd$Orchard, nd$Plant))]
f <- design_anova(nd, "y", c("Orchard", "Plant"), design = "nested")
a <- summary(aov(y ~ Orchard + Error(Orchard:Plant), nd))
chk("Nested: df B(A) = a(b-1)", row_of(f, "Orchard:Plant", "df"), 8)
chk("Nested: F A vs aov Error(A:B)", row_of(f, "Orchard", "F"), a[["Error: Orchard:Plant"]][[1]]["Orchard", "F value"])

# --- Incomplete blocks: alpha-lattice built by agricolae, intra-block analysis vs lm
set.seed(7)
al <- agricolae::design.alpha(paste0("G", 1:12), k = 3, r = 2, seed = 7)$book
names(al)[names(al) == 'paste0("G", 1:12)'] <- "Genotype"
blk <- as.integer(factor(paste(al$replication, al$block)))
al$y <- rnorm(nrow(al), 5 + as.integer(factor(al$Genotype)) / 4) + rnorm(max(blk))[blk]
f <- design_anova(al, "y", "Genotype", block = c("replication", "block"))
al2 <- fac(al, c("replication", "block", "Genotype"))
red <- lm(y ~ replication + replication:block, al2); full <- lm(y ~ replication + replication:block + Genotype, al2)
chk("IBD: suggested design", as.numeric(f$design == "ibd"), 1)
chk("IBD: adjusted SS genotype vs lm", row_of(f, "Genotype", "ss"), anova(red, full)[2, "Sum of Sq"])
chk("IBD: df error vs lm", f$residual[["df"]], df.residual(full))

# --- ANCOVA: maize with a covariate
cv <- maize_nitrogen; set.seed(3); cv$Stand <- round(rnorm(16, 60, 4))
f <- design_anova(cv, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block", covariate = "Stand")
d <- fac(cv, c("Block", "Nitrogen_kg_ha"))
lmc <- lm(Yield_t_ha ~ Stand + Block + Nitrogen_kg_ha, d, contrasts = list(Block = contr.sum, Nitrogen_kg_ha = contr.sum))
chk("ANCOVA: SS III vs car::Anova", f$table$ss, car::Anova(lmc, type = 3)[c("Stand", "Block", "Nitrogen_kg_ha"), "Sum Sq"])
em <- summary(emmeans::emmeans(lmc, "Nitrogen_kg_ha"))
chk("ANCOVA: adjusted means vs emmeans", ls_means(f, "Nitrogen_kg_ha")$mean, em$emmean)
chk("ANCOVA: SE adjusted means vs emmeans", ls_means(f, "Nitrogen_kg_ha")$se, em$SE)

cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
