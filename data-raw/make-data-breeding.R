# Builds the simulated example data of the plant-breeding analyses
# (data/maize_met.rda, maize_diallel.rda, maize_lxt.rda, nc_designs.rda).
# Run from the package root:  Rscript data-raw/make-data-breeding.R
# Every value is simulated here; no data are taken from other packages.

# ---- multi-environment trial: 12 hybrids x 6 environments x 3 blocks ----
set.seed(20261022)
g <- 12; e <- 6; r <- 3
hyb <- sprintf("H%02d", seq_len(g)); env <- paste0("E", seq_len(e))
gen_eff <- rnorm(g, 0, 0.55)
env_eff <- c(-1.9, -0.8, 0.1, 0.6, 1.1, 0.9)
u1 <- scale(rnorm(g))[, 1]; v1 <- scale(rnorm(e))[, 1]
u2 <- scale(rnorm(g))[, 1]; v2 <- scale(rnorm(e))[, 1]
ge <- 0.38 * outer(u1, v1) + 0.2 * outer(u2, v2) + matrix(rnorm(g * e, 0, 0.12), g, e)
maize_met <- expand.grid(Hybrid = hyb, Block = seq_len(r), Environment = env, stringsAsFactors = FALSE)
blk <- matrix(rnorm(e * r, 0, 0.15), e, r)
gi <- match(maize_met$Hybrid, hyb); ei <- match(maize_met$Environment, env)
maize_met$Yield_t_ha <- round(7.2 + gen_eff[gi] + env_eff[ei] + ge[cbind(gi, ei)] +
                                blk[cbind(ei, maize_met$Block)] + rnorm(nrow(maize_met), 0, 0.32), 2)
maize_met <- maize_met[, c("Environment", "Block", "Hybrid", "Yield_t_ha")]

# ---- full diallel among six lines (Griffing method 1), 3 blocks ----
set.seed(77)
p <- 6; par_ <- paste0("P", seq_len(p))
gca <- c(0.45, -0.30, 0.15, -0.55, 0.35, -0.10)
sca <- matrix(0, p, p)
for (i in 1:(p - 1)) for (j in (i + 1):p) sca[i, j] <- sca[j, i] <- rnorm(1, 0.25, 0.35)
mat <- c(0.12, -0.05, 0, -0.10, 0.08, -0.05)
cells <- expand.grid(Female = seq_len(p), Male = seq_len(p))
mu_cell <- with(cells, 6 + gca[Female] + gca[Male] + ifelse(Female == Male, -1.6, 1.2 + sca[cbind(Female, Male)]) +
                  ifelse(Female == Male, 0, mat[Female] - mat[Male] + rnorm(nrow(cells), 0, 0.06)))
maize_diallel <- do.call(rbind, lapply(1:3, function(b) data.frame(Block = b, Female = par_[cells$Female],
                                                                    Male = par_[cells$Male], stringsAsFactors = FALSE)))
maize_diallel$Yield_t_ha <- round(rep(mu_cell, 3) + rep(c(-0.15, 0.05, 0.1), each = nrow(cells)) +
                                    rnorm(nrow(maize_diallel), 0, 0.28), 2)

# ---- line x tester: 6 lines x 3 testers + 9 parents, 3 blocks ----
set.seed(311)
L <- paste0("L", 1:6); Tt <- paste0("T", 1:3)
gl <- c(0.6, -0.4, 0.2, -0.7, 0.4, -0.1); gt <- c(0.35, -0.25, -0.1)
S <- matrix(rnorm(18, 0, 0.25), 6, 3)
ent <- rbind(data.frame(Line = rep(L, 3), Tester = rep(Tt, each = 6), stringsAsFactors = FALSE),
             data.frame(Line = L, Tester = NA_character_, stringsAsFactors = FALSE),
             data.frame(Line = NA_character_, Tester = Tt, stringsAsFactors = FALSE))
mu_ent <- with(ent, ifelse(!is.na(Line) & !is.na(Tester),
                           7 + gl[match(Line, L)] + gt[match(Tester, Tt)] + S[cbind(match(Line, L), match(Tester, Tt))],
                           ifelse(is.na(Tester), 4.6 + 0.6 * gl[match(Line, L)] + rnorm(nrow(ent), 0, 0.2),
                                  4.8 + 0.6 * gt[match(Tester, Tt)] + rnorm(nrow(ent), 0, 0.2))))
maize_lxt <- do.call(rbind, lapply(1:3, function(b) data.frame(Block = b, ent, stringsAsFactors = FALSE)))
maize_lxt$Yield_t_ha <- round(rep(mu_ent, 3) + rep(c(0.1, -0.08, 0), each = nrow(ent)) +
                                rnorm(nrow(maize_lxt), 0, 0.22), 2)

# ---- North Carolina designs, 2 sets x 2 replicates ----
set.seed(4242)
# design I: 3 males per set, 4 females per male
nc1 <- expand.grid(Female = 1:4, Male = 1:3, Rep = 1:2, Set = 1:2)
nc1$Male <- paste0("M", (nc1$Set - 1) * 3 + nc1$Male)
mE <- stats::setNames(rnorm(6, 0, 0.45), paste0("M", 1:6))
fE <- matrix(rnorm(24, 0, 0.5), 6, 4, dimnames = list(paste0("M", 1:6), NULL))
rE <- matrix(rnorm(4, 0, 0.15), 2, 2)
nc1$Yield_t_ha <- round(6.5 + mE[nc1$Male] + fE[cbind(match(nc1$Male, rownames(fE)), nc1$Female)] + rE[cbind(nc1$Set, nc1$Rep)] +
                          rnorm(nrow(nc1), 0, 0.3), 2)
nc1$Female <- paste0("F", nc1$Female)
nc1 <- nc1[, c("Set", "Rep", "Male", "Female", "Yield_t_ha")]
# design II: 4 males x 4 females per set
nc2 <- expand.grid(Female = 1:4, Male = 1:4, Rep = 1:2, Set = 1:2)
nc2$Male <- paste0("M", (nc2$Set - 1) * 4 + nc2$Male); nc2$Female <- paste0("F", (nc2$Set - 1) * 4 + nc2$Female)
m2 <- stats::setNames(rnorm(8, 0, 0.4), paste0("M", 1:8)); f2 <- stats::setNames(rnorm(8, 0, 0.45), paste0("F", 1:8))
mf2 <- matrix(rnorm(64, 0, 0.3), 8, 8, dimnames = list(names(m2), names(f2)))
nc2$Yield_t_ha <- round(6.8 + m2[nc2$Male] + f2[nc2$Female] + mf2[cbind(nc2$Male, nc2$Female)] +
                          rE[cbind(nc2$Set, nc2$Rep)] + rnorm(nrow(nc2), 0, 0.3), 2)
nc2 <- nc2[, c("Set", "Rep", "Male", "Female", "Yield_t_ha")]
# design III: 6 F2 plants per set crossed with both parents; additive a and dominance d at 8 loci
nc3 <- expand.grid(Tester = c("P1", "P2"), Male = 1:6, Rep = 1:2, Set = 1:2, stringsAsFactors = FALSE)
nc3$Male <- paste0("F2-", (nc3$Set - 1) * 6 + nc3$Male)
a <- rep(0.25, 8); dd <- rep(0.18, 8)
geno <- matrix(sample(0:2, 12 * 8, replace = TRUE, prob = c(0.25, 0.5, 0.25)), 12, 8,
               dimnames = list(paste0("F2-", 1:12), NULL))
val <- function(G, tester) {
  # mean of the backcross family: P1 = AA (2), P2 = aa (0); G = number of A alleles of the F2 plant
  pA <- G / 2
  if (tester == "P1") sum(pA * a + (1 - pA) * dd) else sum(pA * dd - (1 - pA) * a)
}
fam <- outer(rownames(geno), c("P1", "P2"), Vectorize(function(i, t) val(geno[i, ], t)))
dimnames(fam) <- list(rownames(geno), c("P1", "P2"))
nc3$Yield_t_ha <- round(6 + fam[cbind(nc3$Male, nc3$Tester)] + rE[cbind(nc3$Set, nc3$Rep)] +
                          rnorm(nrow(nc3), 0, 0.3), 2)
nc3 <- nc3[, c("Set", "Rep", "Male", "Tester", "Yield_t_ha")]
nc_designs <- list(nc1 = nc1, nc2 = nc2, nc3 = nc3)

for (nm in c("maize_met", "maize_diallel", "maize_lxt", "nc_designs"))
  save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz", version = 2)
