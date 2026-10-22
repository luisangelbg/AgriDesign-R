# Stage 5 validation: randomised layouts (structure over many seeds, uniformity,
# reproducibility), field book -> design_anova() round trip against the skeleton
# ANOVA, and replicates / power against power.t.test and the web application.
suppressMessages(devtools::load_all(quiet = TRUE))
out <- "validation/figures"; dir.create(out, showWarnings = FALSE)
ok <- 0; bad <- 0
chk <- function(what, pass, info = "") {
  pass <- isTRUE(all(pass))
  if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %s %s\n", if (pass) "OK" else "FAIL", what, info))
}
all_seeds <- function(n, f) all(vapply(seq_len(n), function(s) isTRUE(all(f(s))), NA))
T4 <- paste0("T", 1:4); T6 <- paste0("G", 1:6)

# ---------------- structure over 200 seeds ----------------
chk("CRD: each treatment r times", all_seeds(200, function(s) {
  b <- make_design("crd", trt = T4, r = 3, seed = s)$book; all(table(b$Treatment) == 3) }))
chk("RCBD: every block holds every treatment once", all_seeds(200, function(s) {
  b <- make_design("rcbd", trt = T6, r = 4, seed = s)$book; all(table(b$Block, b$Treatment) == 1) }))
chk("Latin square: each treatment once per row and per column", all_seeds(200, function(s) {
  b <- make_design("latin", trt = T6, seed = s)$book
  all(table(b$Row, b$Treatment) == 1) && all(table(b$Column, b$Treatment) == 1) && all(b$Row == b$FieldRow) }))
chk("Factorial RCBD: every combination once per block, A and B consistent", all_seeds(200, function(s) {
  b <- make_design("fact_rcbd", A = c("a1", "a2", "a3"), B = c("b1", "b2"), C = c("c1", "c2"), r = 2, seed = s)$book
  all(table(b$Block, b$Treatment) == 1) && all(b$Treatment == paste(b$A, b$B, b$C, sep = " x ")) }))
chk("Split-plot: each main plot one field column holding every B", all_seeds(200, function(s) {
  b <- make_design("split_rcbd", A = c("I1", "I2", "I3"), B = paste0("V", 1:4), r = 3, seed = s)$book
  all(tapply(b$FieldCol, list(b$Block, b$A), function(v) length(unique(v))) == 1) &&
    all(table(paste(b$Block, b$A), b$B) == 1) }))
chk("Strip-plot: A in field rows, B in field columns within each block", all_seeds(200, function(s) {
  b <- make_design("strip", A = c("I1", "I2", "I3"), B = paste0("V", 1:4), r = 3, seed = s)$book
  all(tapply(b$FieldRow, list(b$Block, b$A), function(v) length(unique(v))) == 1) &&
    all(tapply(b$FieldCol, list(b$Block, b$B), function(v) length(unique(v))) == 1) }))
chk("Split-split-plot: every sub-plot one column holding every C", all_seeds(200, function(s) {
  b <- make_design("splitsplit", A = c("a1", "a2"), B = c("b1", "b2", "b3"), C = c("c1", "c2"), r = 2, seed = s)$book
  all(tapply(b$FieldCol, list(paste(b$Block, b$A, b$B)), function(v) length(unique(v))) == 1) &&
    all(table(paste(b$Block, b$A, b$B), b$C) == 1) &&
    all(tapply(b$FieldCol, paste(b$Block, b$A), function(v) diff(range(v))) == 2) }))
chk("Incomplete blocks: every replicate complete, blocks of k plots", all_seeds(200, function(s) {
  b <- make_design("ibd", trt = paste0("G", 1:12), r = 3, k = 4, seed = s)$book
  all(table(b$Block, b$Treatment) == 1) && all(table(paste(b$Block, b$IncBlock)) == 4) }))
chk("Augmented: checks in every block, each entry once", all_seeds(200, function(s) {
  b <- make_design("augmented", checks = c("C1", "C2", "C3"), entries = paste0("E", 1:20), r = 4, seed = s)$book
  all(table(b$Block[b$Type == "check"], b$Treatment[b$Type == "check"]) == 1) &&
    all(table(b$Treatment[b$Type == "entry"]) == 1) && sum(b$Type == "entry") == 20 }))
chk("Plot numbers are 1..n and serpentine neighbours touch", all_seeds(100, function(s) {
  b <- make_design("rcbd", trt = T6, r = 4, seed = s)$book
  b <- b[order(b$Plot), ]
  identical(b$Plot, seq_len(24)) && all(abs(diff(b$FieldRow)) + abs(diff(b$FieldCol)) == 1) }))
b <- make_design("rcbd", trt = T4, r = 3, seed = 5, numbering = "block")$book
chk("Block numbering 101.., 201.., 301..", setequal(b$Plot, c(101:104, 201:204, 301:304)) &&
      all(b$Plot %/% 100 == b$Block))

# ---------------- randomisation: reproducible, uniform, leaves the session RNG alone ----------------
chk("Same seed, same layout", identical(make_design("rcbd", trt = T6, r = 4, seed = 99)$book,
                                        make_design("rcbd", trt = T6, r = 4, seed = 99)$book))
set.seed(123); a1 <- runif(3); set.seed(123); invisible(make_design("rcbd", trt = T6, r = 4, seed = 1)); a2 <- runif(3)
chk("The session random state is restored", identical(a1, a2))
pos <- vapply(1:3000, function(s) { b <- make_design("rcbd", trt = T6, r = 2, seed = s, serpentine = FALSE)$book
                                    b$FieldCol[b$Block == 1 & b$Treatment == "G1"] }, 0)
p <- chisq.test(table(factor(pos, levels = 1:6)))$p.value
chk("RCBD: G1 equally likely in every position (3000 seeds)", p > 0.001, sprintf("(chi-square p = %.3f)", p))
# Latin squares of order 4: there are 576; permuting a cyclic square reaches only 432 of them.
# With the Jacobson-Matthews chain every square must appear, with equal frequency.
sq <- vapply(1:11520, function(s) paste(make_design("latin", trt = T4, seed = s)$book$Treatment, collapse = ""), "")
tb <- table(sq)
chk("Latin 4 x 4: all 576 squares are reached", length(tb) == 576, sprintf("(%d different in 11520 seeds)", length(tb)))
p <- chisq.test(as.vector(tb))$p.value
chk("Latin 4 x 4: squares equally likely", p > 0.001, sprintf("(chi-square p = %.3f, 576 cells, 20 expected each)", p))
is_latin <- function(m) all(apply(m, 1, function(v) length(unique(v)) == length(v))) && all(apply(m, 2, function(v) length(unique(v)) == length(v)))
chk("Jacobson-Matthews squares are Latin for t = 2..12", all(vapply(2:12, function(t) is_latin(AgriDesign:::.with_seed(t, AgriDesign:::.latin_square(t))), NA)))

# ---------------- field book -> design_anova(): df of the skeleton ----------------
round_trip <- function(lay, ...) {
  b <- field_book(lay, "Yield"); set.seed(1); b$Yield <- rnorm(nrow(b), 10)
  f <- design_anova(b, "Yield", ...)
  c(f$table$df[!grepl("^Block$|^Replicates$", f$table$source) | TRUE], f$residual[["df"]])
}
skel <- function(lay) { s <- lay$skeleton; s$df[s$source != "Total"] }
lay <- make_design("crd", trt = T4, r = 4, seed = 1)
chk("Round trip CRD", identical(round_trip(lay, "Treatment"), skel(lay)))
lay <- make_design("rcbd", trt = T6, r = 4, seed = 1)
chk("Round trip RCBD", identical(round_trip(lay, "Treatment", block = "Block"), skel(lay)))
lay <- make_design("latin", trt = T6, seed = 1)
chk("Round trip Latin square", identical(round_trip(lay, "Treatment", row = "Row", col = "Column"), skel(lay)))
lay <- make_design("fact_rcbd", A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 3, seed = 1)
chk("Round trip factorial RCBD", identical(round_trip(lay, c("A", "B"), block = "Block"), skel(lay)))
lay <- make_design("split_rcbd", A = c("I1", "I2"), B = paste0("V", 1:4), r = 3, seed = 1)
chk("Round trip split-plot", identical(round_trip(lay, c("A", "B"), block = "Block", design = "split_rcbd"), skel(lay)))
lay <- make_design("strip", A = c("I1", "I2", "I3"), B = paste0("V", 1:4), r = 3, seed = 1)
rt <- round_trip(lay, c("A", "B"), block = "Block", design = "strip_rcbd")
chk("Round trip strip-plot (same df, other order)", identical(sort(rt), sort(skel(lay))))
lay <- make_design("splitsplit", A = c("a1", "a2"), B = c("b1", "b2", "b3"), C = c("c1", "c2"), r = 3, seed = 1)
chk("Round trip split-split-plot", identical(round_trip(lay, c("A", "B", "C"), block = "Block", design = "splitsplit_rcbd"), skel(lay)))
lay <- make_design("ibd", trt = paste0("G", 1:12), r = 3, k = 4, seed = 1)
chk("Round trip incomplete blocks", identical(round_trip(lay, "Treatment", block = c("Block", "IncBlock")), skel(lay)))
lay <- make_design("augmented", checks = c("C1", "C2", "C3"), entries = paste0("E", 1:12), r = 4, seed = 1)
b <- field_book(lay, "Yield"); set.seed(1); b$Yield <- rnorm(nrow(b), 10)
f <- design_anova(b, "Yield", "Treatment", block = "Block", design = "augmented")
chk("Round trip augmented: error df from the checks", f$residual[["df"]] == lay$skeleton$df[3])

# ---------------- power ----------------
for (r in c(3, 5, 10)) {
  pd <- power_design(cv = 12, diff = 15, t = 2, design = "crd", r = r)$table$power
  ref <- power.t.test(n = r, delta = 15, sd = 12, strict = TRUE)$power
  chk(sprintf("Power, 2 treatments CRD, r = %d, vs power.t.test(strict)", r), abs(pd - ref) < 1e-10, sprintf("(%.6f)", pd))
}
pd <- power_design(cv = 12, diff = 15, t = 2, design = "crd", r = 6)
ref <- power.t.test(n = 6, sd = 12, power = 0.8, strict = TRUE, tol = 1e-10)$delta
chk("Detectable difference vs power.t.test(strict)", abs(pd$table$detectable[pd$table$r == 6] - ref) < 1e-6, sprintf("(%.4f)", ref))
app <- data.frame(cv = c(5, 8, 10, 10, 15, 15, 20), d = c(10, 15, 15, 10, 15, 10, 20), r_app = c(5, 6, 8, 17, 17, 36, 17))
app$r_exact <- mapply(function(cv, d) power_design(cv, d, t = 4)$replicates, app$cv, app$d)
cat("\nReplicates for power 0.80 (4 treatments, RCBD): app (shifted central t) vs exact noncentral t\n")
print(app, row.names = FALSE)
p4 <- power_design(cv = 15, diff = 15, t = 4, r = 4)$table
cat(sprintf("CV 15, d 15, r = 4: power %.3f (app 0.21), detectable %.1f %% (app 33.4 %%)\n", p4$power, p4$detectable))
chk("Exact replicates within one of the app's", all(abs(app$r_exact - app$r_app) <= 1))

# ---------------- figures ----------------
fig <- function(name, expr, w = 16, h = 11) save_figure(file.path(out, paste0(name, ".png")), expr, width = w, height = h, res = 150)
fig("20-layout-rcbd", plot(make_design("rcbd", trt = T6, r = 4, seed = 2026)))
fig("21-layout-split", plot(make_design("split_rcbd", A = c("Drip", "Furrow"), B = paste0("V", 1:4), r = 3, seed = 2026, arrangement = "side", blocks_per_row = 3)))
fig("22-layout-strip", plot(make_design("strip", A = c("Drip", "Furrow", "Sprinkler"), B = paste0("V", 1:4), r = 2, seed = 2026, arrangement = "side", gap = TRUE)))
fig("23-layout-latin", plot(make_design("latin", trt = LETTERS[1:5], seed = 2026)))
fig("24-layout-augmented", plot(make_design("augmented", checks = c("Check A", "Check B"), entries = paste0("G", 1:16), r = 4, seed = 2026)))
fig("25-layout-ibd", plot(make_design("ibd", trt = paste0("G", 1:9), r = 2, k = 3, seed = 2026, arrangement = "side")))
fig("26-power", plot(power_design(cv = 15, diff = 15, t = 4)), w = 20, h = 9)
cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
