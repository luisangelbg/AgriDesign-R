# Reproduces the column coding of the web app (js/lm.js, LM.modelMatrix): every
# term is the product of the sum-to-zero contrasts of its factors, whatever the
# marginal terms in the model. Compares it with R's marginality-aware coding for
# the two designs whose error term has no marginal partner: split-plot in CRD
# (A + A:R + B + A:B) and nested (A + B(A)).
suppressMessages(devtools::load_all(quiet = TRUE))
eff <- function(x) { x <- factor(x); m <- contr.sum(nlevels(x)); m[as.integer(x), , drop = FALSE] }
prod_cols <- function(...) {
  L <- list(...); out <- matrix(1, nrow(L[[1]]), 1)
  for (m in L) out <- do.call(cbind, lapply(seq_len(ncol(out)), function(i) out[, i] * m))
  out
}
app_anova <- function(y, blocks) {           # blocks: list of term matrices
  X <- cbind(1, do.call(cbind, blocks)); asg <- c(0, rep(seq_along(blocks), vapply(blocks, ncol, 0)))
  full <- lm.fit(X, y); sse <- sum(full$residuals^2); dfr <- length(y) - full$rank
  t(vapply(seq_along(blocks), function(i) {
    red <- lm.fit(X[, asg != i, drop = FALSE], y)
    c(df = (length(y) - red$rank) - dfr, ss = sum(red$residuals^2) - sse)
  }, c(0, 0))) |> cbind(df_res = dfr, ms_res = sse / dfr)
}

d <- split_irrigation
A <- eff(d$Irrigation); R <- eff(d$Block); B <- eff(d$Variety)
app <- app_anova(d$Yield_t_ha, list(A = A, AR = prod_cols(A, R), B = B, AB = prod_cols(A, B)))
pkg <- design_anova(d, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_crd")$table
cat("Split-plot in CRD, Error a = replicates within Irrigation\n")
cat(sprintf("  app coding: df Error a = %d, F Irrigation = %.3f (residual df %d)\n",
            app[2, "df"], (app[1, "ss"] / app[1, "df"]) / (app[2, "ss"] / app[2, "df"]), app[1, "df_res"]))
cat(sprintf("  R / package: df Error a = %d, F Irrigation = %.3f\n",
            pkg$df[pkg$term == "Irrigation:Block"], pkg$F[pkg$term == "Irrigation"]))

set.seed(1)
nd <- expand.grid(rep = 1:3, Plant = c("p1", "p2", "p3"), Orchard = c("O1", "O2", "O3", "O4"))
nd$y <- rnorm(nrow(nd), 10 + as.integer(nd$Orchard)) + rnorm(12)[as.integer(interaction(nd$Orchard, nd$Plant))]
A <- eff(nd$Orchard); B <- eff(nd$Plant)
app <- app_anova(nd$y, list(A = A, BA = prod_cols(A, B)))
pkg <- design_anova(nd, "y", c("Orchard", "Plant"), design = "nested")$table
cat("Nested, plant labels repeated within each orchard\n")
cat(sprintf("  app coding: df B(A) = %d, F A = %.3f\n", app[2, "df"],
            (app[1, "ss"] / app[1, "df"]) / (app[2, "ss"] / app[2, "df"])))
cat(sprintf("  R / package: df B(A) = %d, F A = %.3f\n",
            pkg$df[pkg$term == "Orchard:Plant"], pkg$F[pkg$term == "Orchard"]))
