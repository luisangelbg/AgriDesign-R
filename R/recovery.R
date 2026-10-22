# Combined intra- and inter-block analysis of incomplete block designs:
# blocks as random effects, variance components by REML, treatment means
# as generalised least-squares estimates.

# Restricted log-likelihood (up to a constant) for gamma = sigma2_block / sigma2,
# with the error variance profiled out.
.reml_profile <- function(gamma, X, Z, y) {
  n <- length(y); p <- ncol(X)
  V <- diag(n) + gamma * tcrossprod(Z)
  L <- chol(V)
  Xs <- backsolve(L, X, transpose = TRUE); ys <- backsolve(L, y, transpose = TRUE)
  q <- qr(Xs)
  if (q$rank < p) stop("the fixed part of the model is not of full rank", call. = FALSE)
  rss <- sum(qr.resid(q, ys)^2)
  logdetV <- 2 * sum(log(diag(L)))
  logdetXVX <- 2 * sum(log(abs(diag(qr.R(q)))))
  -0.5 * (logdetV + logdetXVX + (n - p) * log(rss / (n - p)))
}

#' Recover inter-block information (REML)
#'
#' In incomplete block designs (alpha, lattice, BIBD, cyclic, augmented),
#' the intra-block analysis of [design_anova()] treats blocks as fixed and
#' uses only comparisons within blocks. When blocks are random, block totals
#' also carry information on the treatments. This function fits blocks as a
#' random effect, estimates the block and error variances by restricted
#' maximum likelihood (REML) and returns the treatment means as generalised
#' least-squares estimates (BLUEs), with their standard errors and the
#' efficiency gained over the intra-block analysis and over a complete block
#' design.
#'
#' @param fit An `"ad_anova"` object of an incomplete block design (designs
#'   `"ibd"`, `"bib"`, `"augmented"`, `"rcbd"` with incomplete blocks).
#' @return An object of class `"ad_recovery"` with the variance components
#'   (`sigma2_block`, `sigma2`), the adjusted treatment means and SEs
#'   (`means`), the average standard error of a difference for the combined,
#'   intra-block and complete-block analyses, and the REML log-likelihood.
#' @examples
#' lay <- make_design("alpha", trt = 12, r = 3, k = 3, seed = 4)
#' b <- field_book(lay, "Yield")
#' b$Yield <- 10 + as.integer(factor(b$Treatment)) / 5 + rep(c(-1, 0, 1), length.out = nrow(b))
#' fit <- design_anova(b, "Yield", "Treatment", block = c("Block", "IncBlock"))
#' recover_interblock(fit)
#' @export
recover_interblock <- function(fit) {
  .check_fit(fit)
  r <- fit$roles; d <- fit$data
  if (!length(r$block)) stop("the design has no blocks", call. = FALSE)
  if (length(r$treatment) != 1 || length(r$covariate)) stop("one treatment factor and no covariates are supported", call. = FALSE)
  trt <- r$treatment
  rep_f <- if (length(r$block) == 2) r$block[1] else NULL
  blk <- if (length(r$block) == 2) interaction(d[[r$block[1]]], d[[r$block[2]]], drop = TRUE) else d[[r$block[1]]]
  fx <- stats::as.formula(paste("~", paste(sprintf("`%s`", c(rep_f, trt)), collapse = " + ")))
  X <- stats::model.matrix(fx, d, contrasts.arg = .contrasts_for(c(rep_f, trt)))
  Z <- stats::model.matrix(~ 0 + blk)
  y <- fit$y
  opt <- stats::optimize(function(lg) -.reml_profile(exp(lg), X, Z, y), c(-12, 8), tol = 1e-10)
  gamma <- exp(opt$minimum)
  bound <- FALSE
  if (.reml_profile(0, X, Z, y) >= -opt$objective) { gamma <- 0; bound <- TRUE }
  n <- length(y); p <- ncol(X)
  V <- diag(n) + gamma * tcrossprod(Z)
  Vi <- chol2inv(chol(V))
  XtViX <- crossprod(X, Vi %*% X)
  cov_b <- solve(XtViX)
  beta <- drop(cov_b %*% crossprod(X, Vi %*% y))
  res <- y - drop(X %*% beta)
  sigma2 <- drop(crossprod(res, Vi %*% res)) / (n - p)
  # treatment means at equal weights of the replicates (as ls_means)
  lv <- levels(d[[trt]])
  grid <- expand.grid(c(if (!is.null(rep_f)) list(levels(d[[rep_f]])), list(lv)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  names(grid) <- c(rep_f, trt)
  for (f in names(grid)) grid[[f]] <- factor(grid[[f]], levels = levels(d[[f]]))
  Xg <- stats::model.matrix(fx, grid, contrasts.arg = .contrasts_for(c(rep_f, trt)))
  L <- t(vapply(lv, function(l) colMeans(Xg[grid[[trt]] == l, , drop = FALSE]), numeric(p)))
  mean_hat <- drop(L %*% beta)
  cov_m <- sigma2 * L %*% cov_b %*% t(L)
  ij <- .pair_index(length(lv))
  sed <- sqrt(cov_m[cbind(ij[, 1], ij[, 1])] + cov_m[cbind(ij[, 2], ij[, 2])] - 2 * cov_m[ij])
  # intra-block analysis and complete-block equivalent for comparison
  intra <- .lsm(fit, trt)
  mse <- fit$residual[["ms"]]
  sed_intra <- vapply(seq_len(nrow(ij)), function(q) {
    dx <- intra$xbar[ij[q, 1], ] - intra$xbar[ij[q, 2], ]; sqrt(drop(t(dx) %*% intra$XtXi %*% dx) * mse) }, 0)
  # complete blocks: the same data analysed with the replicates as complete blocks
  reps <- mean(table(d[[trt]]))
  sed_rcbd <- if (!is.null(rep_f)) {
    f0 <- .ls(stats::model.matrix(fx, d, contrasts.arg = .contrasts_for(c(rep_f, trt))), y)
    sqrt(2 * f0$sse / f0$df / reps)
  } else NA_real_
  out <- data.frame(level = lv, n = as.vector(table(d[[trt]])[lv]), mean = mean_hat, se = sqrt(diag(cov_m)),
                    intra_mean = intra$means$mean, stringsAsFactors = FALSE)
  if (!is.null(fit$transform)) out$back_transformed <- fit$transform$inv(out$mean)
  structure(list(means = out, sigma2_block = gamma * sigma2, sigma2 = sigma2, gamma = gamma, boundary = bound,
                 loglik_reml = -opt$objective, sed = c(combined = mean(sed), intra_block = mean(sed_intra),
                                                       complete_blocks = sed_rcbd),
                 response = fit$response, factor = trt, design_name = fit$design_name),
            class = "ad_recovery")
}

#' @export
print.ad_recovery <- function(x, digits = 4, ...) {
  cat("\nRecovery of inter-block information (REML)\n", x$design_name, "\n", sep = "")
  cat(sprintf("Variance components: blocks = %s, error = %s%s\n", format(x$sigma2_block, digits = digits),
              format(x$sigma2, digits = digits), if (x$boundary) " (block variance on the boundary: no inter-block information)" else ""))
  s <- x$sed
  cat(sprintf("Average SE of a difference: combined %s, intra-block %s%s\n",
              format(s[["combined"]], digits = digits), format(s[["intra_block"]], digits = digits),
              if (is.na(s[["complete_blocks"]])) "" else paste0(", replicates as complete blocks ", format(s[["complete_blocks"]], digits = digits))))
  cat(sprintf("Relative efficiency of the combined analysis: %.1f %% of the intra-block analysis%s\n",
              100 * (s[["intra_block"]] / s[["combined"]])^2,
              if (is.na(s[["complete_blocks"]])) "" else sprintf(", %.1f %% of complete blocks", 100 * (s[["complete_blocks"]] / s[["combined"]])^2)))
  cat("\n"); print(x$means, row.names = FALSE, digits = digits)
  invisible(x)
}
