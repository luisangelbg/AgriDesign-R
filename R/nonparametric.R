# Rank-based and unequal-variance alternatives to the analysis of variance,
# with their pairwise comparisons and compact letter display.

# ---------- compact letter display ----------
# Levels that do not differ are joined in a graph; every maximal clique
# (Bron-Kerbosch) gets a letter. Levels are ordered by `center`, highest
# first, so "a" is the top group; a level may carry several letters when the
# results are not transitive.
.letters <- function(center, different) {
  k <- length(center)
  if (k == 1) return("a")
  ord <- order(-center)
  adj <- !different; diag(adj) <- FALSE
  cliques <- list()
  bk <- function(R, P, X) {
    if (!length(P) && !length(X)) { cliques[[length(cliques) + 1]] <<- R; return(invisible()) }
    for (v in P) {
      bk(c(R, v), P[adj[v, P]], X[adj[v, X]])
      P <- setdiff(P, v); X <- c(X, v)
    }
  }
  bk(integer(0), ord, integer(0))
  rank <- integer(k); rank[ord] <- seq_len(k)
  best <- vapply(cliques, function(cl) min(rank[cl]), 0)
  size <- vapply(cliques, length, 0)
  cliques <- cliques[order(best, -size)]
  lab <- character(k)
  for (ci in seq_along(cliques)) {
    L <- if (ci <= 26) letters[ci] else paste0(letters[(ci - 1) %/% 26], letters[(ci - 1) %% 26 + 1])
    lab[cliques[[ci]]] <- paste0(lab[cliques[[ci]]], L)
  }
  vapply(strsplit(lab, ""), function(s) paste(sort(s), collapse = ""), "")
}

.pair_index <- function(k) {
  ij <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)
  ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
}

.diff_matrix <- function(k, ij, padj, alpha) {
  m <- matrix(FALSE, k, k)
  sig <- padj < alpha
  m[ij] <- sig; m[ij[, 2:1, drop = FALSE]] <- sig
  m
}

.tie_term <- function(x) { t <- table(x); sum(t^3 - t) }

# ---------- tests ----------

.kruskal <- function(y, g) {
  N <- length(y); k <- nlevels(g)
  r <- rank(y)
  mr <- tapply(r, g, mean); n <- tabulate(g)
  H <- 12 / (N * (N + 1)) * sum(n * (mr - (N + 1) / 2)^2)
  C <- 1 - .tie_term(y) / (N^3 - N)
  if (C > 0) H <- H / C
  list(H = H, df = k - 1, p = stats::pchisq(H, k - 1, lower.tail = FALSE),
       epsilon2 = max(0, (H - k + 1) / (N - k)), mean_ranks = mr, N = N, n = n, ties = .tie_term(y))
}

.dunn <- function(kw, ij) {
  N <- kw$N; n <- kw$n
  se <- sqrt((N * (N + 1) / 12 - kw$ties / (12 * (N - 1))) * (1 / n[ij[, 1]] + 1 / n[ij[, 2]]))
  z <- (kw$mean_ranks[ij[, 1]] - kw$mean_ranks[ij[, 2]]) / se
  list(stat = z, p = 2 * stats::pnorm(-abs(z)))
}

# Mann-Whitney by pairs, normal approximation with tie correction and no continuity correction.
.mann_whitney <- function(a, b) {
  n1 <- length(a); n2 <- length(b); N <- n1 + n2
  r <- rank(c(a, b))
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  s <- sqrt(n1 * n2 / 12 * ((N + 1) - .tie_term(c(a, b)) / (N * (N - 1))))
  z <- (U - n1 * n2 / 2) / s
  c(z = z, p = 2 * stats::pnorm(-abs(z)))
}

.games_howell <- function(y, g, ij) {
  n <- tabulate(g); m <- tapply(y, g, mean); v <- tapply(y, g, stats::var)
  a <- v[ij[, 1]] / n[ij[, 1]]; b <- v[ij[, 2]] / n[ij[, 2]]
  se <- sqrt(a + b)
  df <- se^4 / (a^2 / (n[ij[, 1]] - 1) + b^2 / (n[ij[, 2]] - 1))
  q <- abs(m[ij[, 1]] - m[ij[, 2]]) / se * sqrt(2)
  p <- stats::ptukey(q, nlevels(g), df, lower.tail = FALSE)
  list(stat = unname(q), p = unname(p), df = unname(df))
}

# Friedman with the tie-corrected statistic (equal to the classic one without ties).
# `tab`: blocks in rows, treatments in columns.
.friedman <- function(tab) {
  b <- nrow(tab); k <- ncol(tab)
  rk <- t(apply(tab, 1, rank))
  R <- colSums(rk); A <- sum(rk^2); B <- sum(R^2) / b
  Cf <- b * k * (k + 1)^2 / 4
  Fr <- (k - 1) * b * (B - Cf) / (A - Cf)
  Ff <- (b - 1) * Fr / (b * (k - 1) - Fr)
  list(chisq = Fr, df = k - 1, p = stats::pchisq(Fr, k - 1, lower.tail = FALSE),
       W = Fr / (b * (k - 1)), F = Ff, df1 = k - 1, df2 = (b - 1) * (k - 1),
       p_F = stats::pf(Ff, k - 1, (b - 1) * (k - 1), lower.tail = FALSE),
       rank_sums = R, mean_ranks = R / b, A = A, B = B, b = b, k = k)
}

.conover <- function(fr, ij) {
  se <- sqrt(2 * fr$b * (fr$A - fr$B) / ((fr$b - 1) * (fr$k - 1)))
  df <- (fr$b - 1) * (fr$k - 1)
  t <- abs(fr$rank_sums[ij[, 1]] - fr$rank_sums[ij[, 2]]) / se
  list(stat = unname(t), p = unname(2 * stats::pt(t, df, lower.tail = FALSE)))
}

.nemenyi <- function(fr, ij) {
  se <- sqrt(fr$k * (fr$k + 1) / (6 * fr$b))
  q <- abs(fr$mean_ranks[ij[, 1]] - fr$mean_ranks[ij[, 2]]) / se * sqrt(2)
  list(stat = unname(q), p = unname(stats::ptukey(q, fr$k, Inf, lower.tail = FALSE)))
}

# Sum of squares and df of term i (Type III, drop-one) and of the residual, for any response.
.type3_term <- function(X, y, assign, i) {
  full <- .ls(X, y); red <- .ls(X, y, which(assign != i))
  c(ss = red$sse - full$sse, df = red$df - full$df, sse = full$sse, dfe = full$df)
}

#' Rank-based and unequal-variance alternatives to the ANOVA
#'
#' When the assumptions fail, compares the treatments without assuming
#' normality (rank-based tests) or without assuming equal variances
#' (Welch). The method follows the design of `fit`:
#'
#' * `"kruskal"`: Kruskal-Wallis with Dunn's pairwise tests (one factor,
#'   no blocks); `"kruskal_mw"` uses pairwise Mann-Whitney tests instead.
#' * `"welch"`: Welch's one-way ANOVA with Games-Howell pairwise tests (one
#'   factor, no blocks; normal data with unequal variances).
#' * `"friedman"`: Friedman's test with Conover's or Nemenyi's pairwise
#'   tests (one factor in complete blocks, one plot per block and treatment).
#' * `"srh"`: Scheirer-Ray-Hare extension of Kruskal-Wallis (two factors,
#'   no blocks).
#' * `"art"`: aligned rank transform ANOVA (any factorial or blocked design
#'   tested against the residual error), with Tukey comparisons on the
#'   aligned ranks of the first treatment factor.
#'
#' Pairwise p-values of Dunn, Mann-Whitney and Conover are adjusted with
#' `p_adjust`; Games-Howell, Nemenyi and Tukey already control the
#' family-wise error. Levels that share a letter do not differ.
#'
#' @inheritParams check_assumptions
#' @param method One of the methods above; `NULL` picks the first one that
#'   fits the design (Kruskal-Wallis, Friedman, Scheirer-Ray-Hare or ART).
#' @param posthoc For Friedman, `"conover"` (default) or `"nemenyi"`.
#' @param p_adjust Adjustment of the pairwise p-values, as in
#'   [stats::p.adjust()].
#' @return An object of class `"ad_nonparametric"` with the test (`test`),
#'   the table of terms for SRH and ART (`table`), the pairwise comparisons
#'   (`pairs`) and the groups with their letters (`groups`).
#' @examples
#' fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
#' nonparametric(fit)
#'
#' fit <- design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block")
#' nonparametric(fit, posthoc = "nemenyi")
#' @export
nonparametric <- function(fit, method = NULL, posthoc = c("conover", "nemenyi"),
                          p_adjust = "holm", alpha = 0.05) {
  if (!inherits(fit, "ad_anova")) stop("'fit' must come from design_anova()", call. = FALSE)
  posthoc <- match.arg(posthoc)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  r <- fit$roles; d <- fit$data; y <- fit$y
  one <- length(r$treatment) == 1
  no_blk <- !length(r$block) && is.null(r$row)
  f1 <- r$treatment[1]; g <- d[[f1]]; lv <- levels(g); k <- length(lv)
  friedman_ok <- one && length(r$block) == 1 && is.null(r$row) && {
    tb <- table(d[[r$block[1]]], g); all(tb == 1)
  }
  avail <- c(if (one && no_blk) c("kruskal", "kruskal_mw", "welch"),
             if (friedman_ok) "friedman",
             if (length(r$treatment) == 2 && no_blk) "srh",
             if (!any(fit$table$is_error) && is.null(r$covariate)) "art")
  if (is.null(method)) {
    if (!length(avail)) stop("no rank-based method fits this design", call. = FALSE)
    method <- avail[1]
  }
  if (!method %in% avail)
    stop("method '", method, "' does not fit this design; available: ",
         if (length(avail)) paste(avail, collapse = ", ") else "none", call. = FALSE)

  ij <- .pair_index(k)
  pairs <- NULL; tab <- NULL; test <- NULL; center_rank <- rep(NA_real_, k); post <- ""
  mk_pairs <- function(stat_name, stat, p, adjust) {
    padj <- if (adjust) stats::p.adjust(p, p_adjust) else p
    data.frame(level_i = lv[ij[, 1]], level_j = lv[ij[, 2]], statistic_name = stat_name,
               statistic = stat, p = p, p_adj = padj, different = padj < alpha,
               stringsAsFactors = FALSE)
  }
  if (method %in% c("kruskal", "kruskal_mw")) {
    kw <- .kruskal(y, g)
    test <- list(name = "Kruskal-Wallis rank sum test", statistic = c(H = kw$H), df = kw$df, p = kw$p,
                 effect = c(epsilon2 = kw$epsilon2))
    center_rank <- kw$mean_ranks
    if (method == "kruskal") {
      dn <- .dunn(kw, ij); pairs <- mk_pairs("z", dn$stat, dn$p, TRUE); post <- "Dunn"
    } else {
      mw <- t(vapply(seq_len(nrow(ij)), function(q) .mann_whitney(y[g == lv[ij[q, 1]]], y[g == lv[ij[q, 2]]]), c(0, 0)))
      pairs <- mk_pairs("z", mw[, 1], mw[, 2], TRUE); post <- "Mann-Whitney"
    }
  } else if (method == "welch") {
    w <- stats::oneway.test(y ~ g, var.equal = FALSE)
    test <- list(name = "Welch one-way ANOVA (unequal variances)", statistic = c(F = unname(w$statistic)),
                 df = unname(w$parameter), p = w$p.value, effect = NULL)
    gh <- .games_howell(y, g, ij); pairs <- mk_pairs("q", gh$stat, gh$p, FALSE); post <- "Games-Howell"
  } else if (method == "friedman") {
    blk <- d[[r$block[1]]]
    tb <- tapply(y, list(blk, g), function(v) v[1])
    fr <- .friedman(tb)
    test <- list(name = "Friedman rank sum test", statistic = c(chisq = fr$chisq), df = fr$df, p = fr$p,
                 effect = c(Kendall_W = fr$W),
                 iman_davenport = c(F = fr$F, df1 = fr$df1, df2 = fr$df2, p = fr$p_F))
    center_rank <- fr$mean_ranks
    if (posthoc == "conover") { cn <- .conover(fr, ij); pairs <- mk_pairs("t", cn$stat, cn$p, TRUE); post <- "Conover" }
    else { nm <- .nemenyi(fr, ij); pairs <- mk_pairs("q", nm$stat, nm$p, FALSE); post <- "Nemenyi" }
  } else if (method == "srh") {
    rk <- rank(y); ms_tot <- stats::var(rk)
    assign <- attr(fit$X, "assign")
    sse_prev <- .ls(fit$X, rk, which(assign == 0))
    rows <- lapply(seq_along(fit$terms), function(i) {
      cur <- .ls(fit$X, rk, which(assign <= i))
      out <- c(df = sse_prev$df - cur$df, ss = sse_prev$sse - cur$sse)
      sse_prev <<- cur
      out
    })
    rows <- do.call(rbind, rows)
    H <- rows[, "ss"] / ms_tot
    tab <- data.frame(term = fit$table$source, df = rows[, "df"], ss = rows[, "ss"], H = H,
                      p = stats::pchisq(H, rows[, "df"], lower.tail = FALSE), stringsAsFactors = FALSE)
    test <- list(name = "Scheirer-Ray-Hare test (H = SS of ranks / total MS of ranks)", statistic = NULL,
                 df = NULL, p = NULL, effect = NULL)
  } else if (method == "art") {
    X <- fit$X; assign <- attr(X, "assign")
    full <- .ls(X, y); fitted_full <- y - full$fit$residuals
    rows <- lapply(seq_along(fit$terms), function(i) {
      red <- .ls(X, y, which(assign != i))
      aligned <- full$fit$residuals + (fitted_full - (y - red$fit$residuals))
      rk <- rank(aligned)
      a <- .type3_term(X, rk, assign, i)
      F <- (a[["ss"]] / a[["df"]]) / (a[["sse"]] / a[["dfe"]])
      list(row = data.frame(term = fit$table$source[i], df = a[["df"]], df_error = a[["dfe"]], F = F,
                            p = stats::pf(F, a[["df"]], a[["dfe"]], lower.tail = FALSE),
                            aligned_sum = sum(aligned), stringsAsFactors = FALSE),
           ranks = rk, mse = a[["sse"]] / a[["dfe"]], dfe = a[["dfe"]])
    })
    tab <- do.call(rbind, lapply(rows, `[[`, "row"))
    test <- list(name = "Aligned rank transform ANOVA (each term aligned and ranked on its own)",
                 statistic = NULL, df = NULL, p = NULL, effect = NULL)
    ti <- which(fit$table$term == f1)
    if (length(ti) && k > 2 && tab$p[ti] < alpha) {
      rk <- rows[[ti]]$ranks; mse <- rows[[ti]]$mse; dfe <- rows[[ti]]$dfe
      m <- tapply(rk, g, mean); n <- tabulate(g)
      q <- abs(m[ij[, 1]] - m[ij[, 2]]) / sqrt(mse / 2 * (1 / n[ij[, 1]] + 1 / n[ij[, 2]]))
      pairs <- mk_pairs("q", unname(q), unname(stats::ptukey(q, k, dfe, lower.tail = FALSE)), FALSE)
      center_rank <- m; post <- "Tukey on aligned ranks"
    }
  }
  groups <- data.frame(level = lv, n = tabulate(g), median = as.numeric(tapply(y, g, stats::median)),
                       mean = as.numeric(tapply(y, g, mean)), mean_rank = as.numeric(center_rank),
                       stringsAsFactors = FALSE)
  if (!is.null(pairs)) {
    groups$group <- .letters(groups$median, .diff_matrix(k, ij, pairs$p_adj, alpha))
    groups <- groups[order(-groups$median), ]
  }
  structure(list(method = method, test = test, table = tab, pairs = pairs, groups = groups,
                 posthoc = post, p_adjust = if (post %in% c("Dunn", "Mann-Whitney", "Conover")) p_adjust else "none",
                 alpha = alpha, factor = f1, response = fit$response, available = avail),
            class = "ad_nonparametric")
}

#' @export
print.ad_nonparametric <- function(x, digits = 4, ...) {
  cat("\n", x$test$name, "\n", sep = "")
  cat("Response: ", x$response, "   factor: ", x$factor, "\n", sep = "")
  if (!is.null(x$test$statistic)) {
    cat(sprintf("%s = %s, df = %s, p = %s\n", names(x$test$statistic), format(x$test$statistic, digits = digits),
                paste(format(x$test$df, digits = digits), collapse = ", "), format.pval(x$test$p, digits = digits)))
    if (!is.null(x$test$effect))
      cat(sprintf("%s = %s\n", names(x$test$effect), format(x$test$effect, digits = digits)))
    if (!is.null(x$test$iman_davenport)) {
      id <- x$test$iman_davenport
      cat(sprintf("Iman-Davenport F = %s on %g and %g df, p = %s\n", format(id[["F"]], digits = digits),
                  id[["df1"]], id[["df2"]], format.pval(id[["p"]], digits = digits)))
    }
  }
  if (!is.null(x$table)) {
    cat("\n"); print(x$table, row.names = FALSE, digits = digits)
  }
  if (!is.null(x$pairs)) {
    cat("\nPairwise comparisons: ", x$posthoc,
        if (x$p_adjust != "none") paste0(" (p adjusted by ", x$p_adjust, ")"), "\n", sep = "")
    print(x$pairs[c("level_i", "level_j", "statistic", "p", "p_adj", "different")], row.names = FALSE, digits = digits)
    cat("\nGroups (levels sharing a letter do not differ at alpha = ", x$alpha, "):\n", sep = "")
    print(x$groups, row.names = FALSE, digits = digits)
  }
  invisible(x)
}
