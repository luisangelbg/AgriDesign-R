# Builds the simulated example data of R/pathology.R and R/tests-extra.R:
#   blight_progress   disease progress of late blight under four fungicide programmes
#   uniformity_trial  uniformity trial of 24 x 16 basic units (Smith's index)
#   bib_sorghum       seven sorghum hybrids in a balanced incomplete block design
# Every value is simulated here; nothing is taken from other packages.
# Run from the package root:  Rscript data-raw/make-data-pathology.R

# ---------- blight_progress ----------
# Logistic progress y = 100 / (1 + exp(-r (t - t50))) with a rate and a
# half-time per programme, a block effect on the half-time and assessment noise.
set.seed(20261022)
trt <- data.frame(Treatment = c("Untreated", "Contact", "Systemic", "Alternated"),
                  r = c(0.21, 0.17, 0.13, 0.12), t50 = c(30, 38, 47, 50))
days <- c(21, 28, 35, 42, 49, 56, 63)
blk_shift <- c(-2.5, 0.8, 1.6, 0.1)
rows <- list()
for (b in 1:4) for (i in seq_len(nrow(trt))) {
  t50 <- trt$t50[i] + blk_shift[b] + stats::rnorm(1, 0, 1.2)
  r <- trt$r[i] * exp(stats::rnorm(1, 0, 0.05))
  y <- 100 / (1 + exp(-r * (days - t50)))
  y <- y * exp(stats::rnorm(length(days), 0, 0.08))
  y <- pmin(100, pmax(0, round(cummax(y), 1)))
  rows[[length(rows) + 1]] <- data.frame(Block = b, Treatment = trt$Treatment[i], Day = days,
                                         Severity_pct = y, stringsAsFactors = FALSE)
}
blight_progress <- do.call(rbind, rows)
rownames(blight_progress) <- NULL

# ---------- uniformity_trial ----------
# Yield (kg) of 24 rows x 16 columns of basic units: a fertility gradient along
# the rows, a smooth patchy component (moving average of white noise) and
# independent noise.
set.seed(1938)
nr <- 24; nc <- 16
grad <- outer(seq(-1, 1, length.out = nr), rep(1, nc)) * 0.35
z <- matrix(stats::rnorm((nr + 4) * (nc + 4)), nr + 4, nc + 4)
patch <- matrix(0, nr, nc)
for (i in 1:nr) for (j in 1:nc) patch[i, j] <- mean(z[i:(i + 4), j:(j + 4)])
uniformity_trial <- round(5 + grad + 0.9 * patch + stats::rnorm(nr * nc, 0, 0.32), 2)
dimnames(uniformity_trial) <- list(paste0("R", 1:nr), paste0("C", 1:nc))

# ---------- bib_sorghum ----------
# BIBD t = 7, k = 3, b = 7, r = 3, lambda = 1 (cyclic from the difference set {0, 1, 3}).
set.seed(7)
hyb <- paste0("H", 1:7)
eff <- c(0.0, 0.45, -0.30, 0.95, 0.20, -0.65, 0.60)
rows <- list()
for (b in 1:7) {
  plots <- ((b - 1 + c(0, 1, 3)) %% 7) + 1
  be <- stats::rnorm(1, 0, 0.5)
  rows[[b]] <- data.frame(Block = b, Hybrid = hyb[plots],
                          Yield_t_ha = round(6 + eff[plots] + be + stats::rnorm(3, 0, 0.3), 2),
                          stringsAsFactors = FALSE)
}
bib_sorghum <- do.call(rbind, rows)
rownames(bib_sorghum) <- NULL

for (nm in c("blight_progress", "uniformity_trial", "bib_sorghum"))
  save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz", version = 2)
