# Randomised field layouts of the classical designs, field books, skeleton
# ANOVA tables, and replicates / power for a difference between two means.

# Runs `expr` with the given seed and restores the random state of the session.
.with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

.shuffle <- function(x) x[sample.int(length(x))]

# Uniformly distributed random Latin square of order t (Jacobson and Matthews,
# 1996): a Markov chain on the incidence cube M[row, col, symbol] that visits
# every Latin square with equal probability. Permuting the rows, columns and
# symbols of a cyclic square only reaches its isotopy class (432 of the 576
# squares of order 4), so it is not a uniform randomisation.
# The square is the state after a fixed number of steps; when that state is
# improper the chain starts again. Waiting for the first proper state after
# the fixed number of steps would bias the draw (checked with 11520 squares
# of order 4: chi-square p < 1e-4 against p = 0.04 with the fixed number).
.latin_square <- function(t, steps = max(1000, 4 * t^3)) {
  repeat {
    M <- .jm_chain(t, steps)
    if (!is.null(M)) return(apply(M, c(1, 2), which.max))
  }
}

# `steps` moves of the chain from a cyclic square; NULL when it ends improper.
.jm_chain <- function(t, steps) {
  M <- array(0L, c(t, t, t))
  for (i in seq_len(t)) for (j in seq_len(t)) M[i, j, (i + j - 2) %% t + 1] <- 1L
  improper <- NULL
  pick <- function(v) if (length(v) == 1) v else sample(v, 1)
  for (n in seq_len(steps)) {
    if (is.null(improper)) {
      repeat { r <- sample.int(t, 1); c <- sample.int(t, 1); s <- sample.int(t, 1); if (M[r, c, s] == 0L) break }
      r2 <- which(M[, c, s] == 1L); c2 <- which(M[r, , s] == 1L); s2 <- which(M[r, c, ] == 1L)
    } else {
      r <- improper[1]; c <- improper[2]; s <- improper[3]
      r2 <- pick(which(M[, c, s] == 1L)); c2 <- pick(which(M[r, , s] == 1L)); s2 <- pick(which(M[r, c, ] == 1L))
    }
    M[r, c, s] <- M[r, c, s] + 1L; M[r, c2, s2] <- M[r, c2, s2] + 1L
    M[r2, c, s2] <- M[r2, c, s2] + 1L; M[r2, c2, s] <- M[r2, c2, s] + 1L
    M[r, c, s2] <- M[r, c, s2] - 1L; M[r, c2, s] <- M[r, c2, s] - 1L
    M[r2, c, s] <- M[r2, c, s] - 1L; M[r2, c2, s2] <- M[r2, c2, s2] - 1L
    improper <- if (M[r2, c2, s2] < 0L) c(r2, c2, s2) else NULL
  }
  if (is.null(improper)) M else NULL
}

# Cells of one block: a data frame of treatment attributes with their row and
# column inside the block.
.to_grid <- function(cells, cols) {
  i <- seq_len(nrow(cells)) - 1
  cells$r <- i %/% cols; cells$c <- i %% cols
  cells
}

# Put the blocks on the field, stacked or side by side; returns the plots with
# their field row and column (from 1) and the outlines of blocks and main plots.
.arrange <- function(blocks, arrangement = "stacked", blocks_per_row = 2, gap = FALSE) {
  per_row <- if (arrangement == "side") blocks_per_row else 1
  r0 <- 0; c0 <- 0; max_rows <- 0
  plots <- list(); boxes <- list()
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    if (bi > 1 && (bi - 1) %% per_row == 0) { r0 <- r0 + max_rows + gap; c0 <- 0; max_rows <- 0 }
    nr <- max(b$cells$r) + 1; nc <- max(b$cells$c) + 1
    p <- b$cells; p$FieldRow <- r0 + p$r + 1; p$FieldCol <- c0 + p$c + 1; p$block_name <- b$name
    plots[[bi]] <- p
    boxes[[length(boxes) + 1]] <- data.frame(name = b$name, r0 = r0 + 1, c0 = c0 + 1, r1 = r0 + nr, c1 = c0 + nc,
                                             kind = "block", stringsAsFactors = FALSE)
    if (!is.null(b$boxes)) {
      bx <- b$boxes
      boxes[[length(boxes) + 1]] <- data.frame(name = bx$name, r0 = r0 + bx$r0 + 1, c0 = c0 + bx$c0 + 1,
                                               r1 = r0 + bx$r1 + 1, c1 = c0 + bx$c1 + 1, kind = bx$kind,
                                               stringsAsFactors = FALSE)
    }
    c0 <- c0 + nc + gap; max_rows <- max(max_rows, nr)
  }
  all <- Reduce(union, lapply(plots, names))
  plots <- do.call(rbind, lapply(plots, function(p) { for (m in setdiff(all, names(p))) p[[m]] <- NA; p[all] }))
  list(plots = plots, boxes = do.call(rbind, boxes))
}

.cells <- function(trt, ...) data.frame(trt = trt, ..., stringsAsFactors = FALSE)

.gen <- list(
  crd = function(o) {
    lst <- .shuffle(rep(o$trt, each = o$r))
    cols <- if (is.null(o$plots_per_row)) ceiling(sqrt(length(lst))) else o$plots_per_row
    list(list(name = "Field", cells = .to_grid(.cells(lst, rep = NA, block_index = 0, in_block = seq_along(lst) - 1), cols)))
  },
  rcbd = function(o) lapply(seq_len(o$r), function(b) {
    p <- .shuffle(o$trt)
    list(name = paste("Block", b), cells = .to_grid(.cells(p, rep = b, block_index = b - 1, in_block = seq_along(p) - 1),
                                                    if (is.null(o$plots_per_row)) length(p) else o$plots_per_row))
  }),
  latin = function(o) {
    t <- length(o$trt)
    sq <- .latin_square(t); sp <- .shuffle(o$trt)
    g <- expand.grid(j = seq_len(t), i = seq_len(t))
    cells <- .cells(sp[sq[cbind(g$i, g$j)]], Row = g$i, Column = g$j, rep = NA, block_index = 0,
                    in_block = (g$i - 1) * t + g$j - 1)
    cells$r <- g$i - 1; cells$c <- g$j - 1
    list(list(name = "Square", cells = cells))
  },
  split = function(o) lapply(seq_len(o$r), function(b) {
    mains <- .shuffle(o$A); nb <- length(o$B)
    cells <- do.call(rbind, lapply(seq_along(mains), function(mi) {
      s <- .shuffle(o$B)
      .cells(paste(mains[mi], s, sep = " x "), A = mains[mi], B = s, rep = b, block_index = b - 1,
             r = seq_len(nb) - 1, c = mi - 1)
    }))
    cells$in_block <- seq_len(nrow(cells)) - 1
    list(name = paste("Block", b), cells = cells,
         boxes = data.frame(name = mains, r0 = 0, c0 = seq_along(mains) - 1, r1 = nb - 1, c1 = seq_along(mains) - 1,
                            kind = "main", stringsAsFactors = FALSE))
  }),
  strip = function(o) lapply(seq_len(o$r), function(b) {
    ra <- .shuffle(o$A); cb <- .shuffle(o$B)
    g <- expand.grid(j = seq_along(cb), i = seq_along(ra))
    cells <- .cells(paste(ra[g$i], cb[g$j], sep = " x "), A = ra[g$i], B = cb[g$j], rep = b, block_index = b - 1,
                    r = g$i - 1, c = g$j - 1, in_block = seq_len(nrow(g)) - 1)
    list(name = paste("Block", b), cells = cells,
         boxes = rbind(data.frame(name = ra, r0 = seq_along(ra) - 1, c0 = 0, r1 = seq_along(ra) - 1, c1 = length(cb) - 1,
                                  kind = "main", stringsAsFactors = FALSE),
                       data.frame(name = cb, r0 = 0, c0 = seq_along(cb) - 1, r1 = length(ra) - 1, c1 = seq_along(cb) - 1,
                                  kind = "sub", stringsAsFactors = FALSE)))
  }),
  splitsplit = function(o) lapply(seq_len(o$r), function(b) {
    mains <- .shuffle(o$A); col <- 0; cells <- list(); bx <- list()
    for (a in mains) {
      c0 <- col
      for (s in .shuffle(o$B)) {
        ss <- .shuffle(o$C)
        cells[[length(cells) + 1]] <- .cells(paste(a, s, ss, sep = " x "), A = a, B = s, C = ss, rep = b,
                                             block_index = b - 1, r = seq_along(ss) - 1, c = col)
        bx[[length(bx) + 1]] <- data.frame(name = s, r0 = 0, c0 = col, r1 = length(ss) - 1, c1 = col, kind = "sub")
        col <- col + 1
      }
      bx[[length(bx) + 1]] <- data.frame(name = a, r0 = 0, c0 = c0, r1 = length(o$C) - 1, c1 = col - 1, kind = "main")
    }
    cells <- do.call(rbind, cells); cells$in_block <- seq_len(nrow(cells)) - 1
    list(name = paste("Block", b), cells = cells, boxes = do.call(rbind, bx))
  }),
  ibd = function(o) {
    t <- length(o$trt); k <- o$k
    if (t %% k) stop("the number of treatments (", t, ") must be a multiple of the block size (", k, ")", call. = FALSE)
    nb <- t / k
    lapply(seq_len(o$r), function(rr) {
      p <- .shuffle(o$trt)
      cells <- .cells(p, rep = rr, IncBlock = rep(seq_len(nb), each = k), block_index = rr - 1,
                      r = rep(seq_len(nb) - 1, each = k), c = rep(seq_len(k) - 1, nb), in_block = seq_len(t) - 1)
      list(name = paste("Rep", rr), cells = cells,
           boxes = data.frame(name = paste0("B", seq_len(nb)), r0 = seq_len(nb) - 1, c0 = 0, r1 = seq_len(nb) - 1,
                              c1 = k - 1, kind = "sub", stringsAsFactors = FALSE))
    })
  },
  augmented = function(o) {
    e <- .shuffle(o$entries); per <- ceiling(length(e) / o$r)
    lapply(seq_len(o$r), function(b) {
      mine <- e[seq_len(per) + (b - 1) * per]; mine <- mine[!is.na(mine)]
      lst <- .shuffle(c(o$checks, mine))
      list(name = paste("Block", b),
           cells = .to_grid(.cells(lst, rep = b, Type = ifelse(lst %in% o$checks, "check", "entry"), block_index = b - 1,
                                   in_block = seq_along(lst) - 1),
                            if (is.null(o$plots_per_row)) length(lst) else o$plots_per_row))
    })
  })

# Skeleton ANOVA (sources and degrees of freedom) of each generated design.
.skeleton <- function(design, t, r, a, b, cc, k, nc, nb = 0) {
  rows <- switch(design,
    graeco = list(c("Rows", t - 1), c("Columns", t - 1), c("Greek letters", t - 1), c("Treatments", t - 1),
                  c("Error", (t - 1) * (t - 3))),
    youden = list(c("Rows", k - 1), c("Columns (blocks)", t - 1), c("Treatments (adjusted)", t - 1),
                  c("Error", t * k - 1 - (k - 1) - 2 * (t - 1))),
    bib = , cyclic = list(c("Blocks", nb - 1), c("Treatments (adjusted)", t - 1), c("Intra-block error", nb * k - nb - t + 1)),
    alpha = , lattice = list(c("Replicates", r - 1), c("Blocks within replicates", r * (t / k - 1)),
                             c("Treatments (adjusted)", t - 1), c("Intra-block error", (r - 1) * (t - 1) - r * (t / k - 1))),
    crd = list(c("Treatments", t - 1), c("Error", t * (r - 1))),
    rcbd = list(c("Blocks", r - 1), c("Treatments", t - 1), c("Error", (r - 1) * (t - 1))),
    latin = list(c("Rows", t - 1), c("Columns", t - 1), c("Treatments", t - 1), c("Error", (t - 1) * (t - 2))),
    fact_crd = , fact_rcbd = {
      ab <- a * b * max(cc, 1)
      x <- c(if (design == "fact_rcbd") list(c("Blocks", r - 1)),
             list(c("A", a - 1), c("B", b - 1), c("A x B", (a - 1) * (b - 1))),
             if (cc) list(c("C", cc - 1), c("A x C", (a - 1) * (cc - 1)), c("B x C", (b - 1) * (cc - 1)),
                          c("A x B x C", (a - 1) * (b - 1) * (cc - 1))))
      c(x, list(c("Error", if (design == "fact_rcbd") (ab - 1) * (r - 1) else ab * (r - 1))))
    },
    split_rcbd = list(c("Blocks", r - 1), c("A (main plot)", a - 1), c("Error a", (r - 1) * (a - 1)),
                      c("B (sub-plot)", b - 1), c("A x B", (a - 1) * (b - 1)), c("Error b", a * (r - 1) * (b - 1))),
    strip = list(c("Blocks", r - 1), c("A", a - 1), c("Error a", (r - 1) * (a - 1)), c("B", b - 1),
                 c("Error b", (r - 1) * (b - 1)), c("A x B", (a - 1) * (b - 1)), c("Error c", (r - 1) * (a - 1) * (b - 1))),
    splitsplit = list(c("Blocks", r - 1), c("A", a - 1), c("Error a", (r - 1) * (a - 1)), c("B", b - 1),
                      c("A x B", (a - 1) * (b - 1)), c("Error b", a * (r - 1) * (b - 1)), c("C", cc - 1),
                      c("A x C", (a - 1) * (cc - 1)), c("B x C", (b - 1) * (cc - 1)),
                      c("A x B x C", (a - 1) * (b - 1) * (cc - 1)), c("Error c", a * b * (r - 1) * (cc - 1))),
    ibd = list(c("Replicates", r - 1), c("Blocks within replicates", r * (t / k - 1)),
               c("Treatments (adjusted)", t - 1), c("Intra-block error", (r - 1) * (t - 1) - r * (t / k - 1))),
    augmented = list(c("Blocks", r - 1), c("Checks", nc - 1), c("Error (from checks)", (r - 1) * (nc - 1))))
  out <- data.frame(source = vapply(rows, `[`, "", 1), df = as.numeric(vapply(rows, `[`, "", 2)), stringsAsFactors = FALSE)
  n <- switch(design, latin = , graeco = t * t, youden = t * k, bib = , cyclic = nb * k,
              fact_crd = , fact_rcbd = a * b * max(cc, 1) * r, split_rcbd = , strip = a * b * r,
              splitsplit = a * b * cc * r, t * r)
  if (design != "augmented") out <- rbind(out, data.frame(source = "Total", df = n - 1))
  out
}

.generated <- c(crd = "Completely randomised design (CRD)", rcbd = "Randomised complete block design (RCBD)",
                latin = "Latin square", fact_crd = "Factorial in CRD", fact_rcbd = "Factorial in RCBD",
                split_rcbd = "Split-plot in RCBD", strip = "Strip-plot (split-block) in RCBD",
                splitsplit = "Split-split-plot in RCBD",
                ibd = "Resolvable incomplete blocks (random blocks of size k within replicates)",
                augmented = "Augmented design (checks in every block + unreplicated entries)",
                graeco = "Graeco-Latin square", youden = "Youden square",
                bib = "Balanced incomplete block design (BIBD)", cyclic = "Cyclic incomplete block design",
                alpha = "Alpha design (resolvable incomplete blocks)", lattice = "Square lattice")

#' Randomise a field experiment
#'
#' Generates the randomised layout of a classical design, numbers the plots
#' and returns the field book (one row per plot) with the skeleton analysis
#' of variance. The randomisation is reproducible with `seed` and does not
#' change the random state of the session.
#'
#' Designs: `"crd"`, `"rcbd"`, `"latin"` (as many rows and columns as
#' treatments), `"fact_crd"` and `"fact_rcbd"` (all combinations of `A`,
#' `B` and optionally `C`), `"split_rcbd"` (`A` on main plots, `B` on
#' sub-plots), `"strip"` (`A` in horizontal and `B` in vertical strips),
#' `"splitsplit"` (`A`, `B`, `C`), `"ibd"` (each replicate split at random
#' into blocks of `k` plots; the number of treatments must be a multiple of
#' `k`) and `"augmented"` (`checks` in every block and each of the
#' `entries` once).
#'
#' Further designs: `"graeco"` (Graeco-Latin square: rows, columns, a
#' second blocking factor `greek` and the treatments, for an odd number of
#' treatments or 4 or 8), `"youden"` (incomplete Latin square: `k` rows and
#' as many columns as treatments, from a cyclic difference set), `"bib"`
#' (balanced incomplete blocks of size `k`: cyclic when a difference set is
#' known, otherwise all combinations), `"cyclic"` (incomplete blocks
#' developed from initial blocks chosen for the highest efficiency; `r` a
#' multiple of `k`), `"alpha"` (alpha design: resolvable, `t / k` blocks of
#' size `k` per replicate, with the alpha array of highest efficiency found)
#' and `"lattice"` (square lattice: `k^2` treatments in blocks of `k`).
#' [design_efficiency()] reports the efficiency factor of the layout.
#'
#' @param design Design identifier, see Details.
#' @param trt Treatment names, or a single number of treatments (named T1,
#'   T2...).
#' @param r Number of replicates (blocks).
#' @param A,B,C Levels of the factors of factorial and split designs.
#' @param k Block size of incomplete block designs.
#' @param checks,entries Checks and new entries of an augmented design.
#' @param greek Levels of the second blocking factor of a Graeco-Latin square.
#' @param seed Seed of the randomisation; `NULL` draws one and reports it.
#' @param serpentine Number the plots in serpentine order (back and forth
#'   along the field rows).
#' @param numbering `"sequential"` (1, 2, 3...) or `"block"` (101, 102...
#'   in block 1, 201... in block 2).
#' @param arrangement Blocks `"stacked"` one below the other or `"side"` by
#'   side.
#' @param blocks_per_row Blocks per row of the field with `arrangement = "side"`.
#' @param plots_per_row Plots per row within a block (CRD, RCBD, factorials,
#'   augmented); by default a whole block in one row.
#' @param gap Leave an empty row or column between blocks.
#' @return An object of class `"ad_layout"` with the field book (`book`),
#'   the outlines of blocks and main plots (`boxes`), the skeleton ANOVA
#'   (`skeleton`) and the seed used.
#' @examples
#' lay <- make_design("rcbd", trt = c("N0", "N60", "N120", "N180"), r = 4, seed = 2026)
#' lay
#' head(lay$book)
#' plot(lay)
#'
#' make_design("split_rcbd", A = c("Drip", "Furrow"), B = paste0("V", 1:4), r = 3, seed = 1)
#' @export
make_design <- function(design = names(.generated), trt = NULL, r = 3, A = NULL, B = NULL, C = NULL, k = NULL,
                        checks = NULL, entries = NULL, greek = NULL, seed = NULL, serpentine = TRUE,
                        numbering = c("sequential", "block"), arrangement = c("stacked", "side"),
                        blocks_per_row = 2, plots_per_row = NULL, gap = FALSE) {
  design <- match.arg(design); numbering <- match.arg(numbering); arrangement <- match.arg(arrangement)
  if (is.numeric(trt) && length(trt) == 1) trt <- paste0("T", seq_len(trt))
  if (!is.null(trt)) trt <- as.character(trt)
  need <- function(x, what, n = 2) if (length(x) < n) stop(design, " needs ", what, call. = FALSE)
  if (design %in% c("crd", "rcbd", "latin", "ibd", "graeco", "youden", "bib", "cyclic", "alpha", "lattice"))
    need(unique(trt), "at least two treatments in 'trt'")
  if (design %in% c("youden", "bib", "cyclic", "alpha", "lattice") && (is.null(k) || k < 2 || k >= length(trt)))
    stop("give a block size 'k' between 2 and the number of treatments - 1", call. = FALSE)
  if (design == "lattice" && length(trt) != k^2) stop("a square lattice needs k^2 treatments", call. = FALSE)
  if (design == "graeco" && is.null(greek)) greek <- paste0("g", seq_along(trt))
  if (design %in% c("fact_crd", "fact_rcbd", "split_rcbd", "strip", "splitsplit")) { need(A, "two or more levels of 'A'"); need(B, "two or more levels of 'B'") }
  if (design == "splitsplit") need(C, "two or more levels of 'C'")
  if (design == "ibd" && (is.null(k) || k < 2 || k >= length(trt))) stop("give a block size 'k' between 2 and the number of treatments - 1", call. = FALSE)
  if (design == "augmented") { need(checks, "two or more 'checks'"); need(entries, "some 'entries'", 1)
    if (length(intersect(checks, entries))) stop("checks and entries must have different names", call. = FALSE) }
  if (any(duplicated(trt)) || any(duplicated(A)) || any(duplicated(B)) || any(duplicated(C)))
    stop("treatment and level names must be unique", call. = FALSE)
  if (!design %in% c("latin", "graeco", "youden", "bib") && r < 2) stop("'r' must be at least 2", call. = FALSE)
  if (is.null(seed)) seed <- sample.int(99999, 1)
  o <- list(trt = trt, r = r, A = A, B = B, C = C, k = k, checks = checks, entries = entries, plots_per_row = plots_per_row,
            greek = as.character(greek), lattice = design == "lattice")
  if (design %in% c("fact_crd", "fact_rcbd")) {
    g <- expand.grid(A = A, B = B, C = if (length(C)) C else NA, stringsAsFactors = FALSE)
    o$trt <- do.call(paste, c(g[if (length(C)) 1:3 else 1:2], sep = " x "))
  }
  gen <- switch(design, fact_crd = .gen$crd, fact_rcbd = .gen$rcbd, split_rcbd = .gen$split, lattice = .gen_extra$alpha,
                c(.gen, .gen_extra)[[design]])
  L <- .with_seed(seed, .arrange(gen(o), arrangement, blocks_per_row, gap))
  p <- L$plots
  if (design %in% c("fact_crd", "fact_rcbd")) {
    m <- match(p$trt, o$trt)
    p$A <- g$A[m]; p$B <- g$B[m]; if (length(C)) p$C <- g$C[m]
  }
  # plot numbers along the field rows (serpentine) or by block
  ord <- order(p$FieldRow, ifelse(serpentine & (match(p$FieldRow, sort(unique(p$FieldRow))) %% 2 == 0), -p$FieldCol, p$FieldCol))
  p$Plot <- NA_integer_
  p$Plot[ord] <- seq_len(nrow(p))
  if (numbering == "block") {
    rank_in <- stats::ave(p$Plot, p$block_index, FUN = rank)
    p$Plot <- as.integer((p$block_index + 1) * 100 + rank_in)
  }
  p <- p[order(p$Plot), ]
  book <- data.frame(Plot = p$Plot, stringsAsFactors = FALSE)
  if (any(!is.na(p$rep))) book$Block <- p$rep
  if ("IncBlock" %in% names(p)) book$IncBlock <- p$IncBlock
  if (design %in% c("latin", "graeco", "youden")) { book$Row <- p$Row; book$Column <- p$Column }
  if (design == "graeco") book$Greek <- p$Greek
  book$FieldRow <- p$FieldRow; book$FieldCol <- p$FieldCol
  if ("A" %in% names(p) && any(!is.na(p$A))) { book$A <- p$A; book$B <- p$B; if (length(C)) book$C <- p$C }
  book$Treatment <- p$trt
  if ("Type" %in% names(p)) book$Type <- p$Type
  rownames(book) <- NULL
  nt <- if (design == "augmented") length(checks) + length(entries) else length(unique(book$Treatment))
  nb <- if ("IncBlock" %in% names(book)) length(unique(paste(book$Block, book$IncBlock)))
        else if ("Block" %in% names(book)) length(unique(book$Block)) else 0
  reps <- if (design %in% c("latin", "graeco", "youden", "bib", "cyclic")) round(nrow(book) / length(unique(book$Treatment))) else r
  skel <- .skeleton(design, length(o$trt), reps, length(A), length(B), length(C), if (is.null(k)) 0 else k, length(checks), nb)
  structure(list(book = book, boxes = L$boxes, skeleton = skel, design = design, design_name = .generated[[design]],
                 seed = seed, r = reps, treatments = nt,
                 size = c(rows = max(book$FieldRow), cols = max(book$FieldCol))),
            class = "ad_layout")
}

#' @export
print.ad_layout <- function(x, ...) {
  cat("\n", x$design_name, "\n", sep = "")
  cat(sprintf("%d plots in %d field rows x %d columns; %d treatments; %s %d; seed %s\n", nrow(x$book), x$size[["rows"]],
              x$size[["cols"]], x$treatments, if (x$design == "latin") "rows = columns =" else "replicates", x$r, x$seed))
  cat("\nSkeleton ANOVA\n"); print(x$skeleton, row.names = FALSE)
  err <- x$skeleton$df[grepl("^(Error|Intra)", x$skeleton$source)]
  if (length(err) && min(err) < 10) cat("Note: fewer than 10 error degrees of freedom; consider more replicates.\n")
  cat("\nField book (first rows):\n"); print(utils::head(x$book, 8), row.names = FALSE)
  invisible(x)
}

#' Field book of a randomised layout
#'
#' The field book of [make_design()] with empty columns for the responses to
#' be measured, ready to be filled in and analysed with [design_anova()].
#'
#' @param layout An `"ad_layout"` object.
#' @param responses Names of the response columns to add.
#' @param file Optional CSV file to write.
#' @return The field book as a data frame (invisibly when `file` is given).
#' @examples
#' lay <- make_design("rcbd", trt = 4, r = 3, seed = 7)
#' head(field_book(lay, c("Yield", "Height")))
#' @export
field_book <- function(layout, responses = "Response", file = NULL) {
  if (!inherits(layout, "ad_layout")) stop("'layout' must come from make_design()", call. = FALSE)
  b <- layout$book
  for (v in responses) b[[v]] <- NA_real_
  if (!is.null(file)) { utils::write.csv(b, file, row.names = FALSE, na = ""); return(invisible(b)) }
  b
}

#' Plot a randomised field layout
#'
#' Draws the plots on their field rows and columns coloured by treatment
#' (or by main-plot factor), with the plot numbers and the outlines of
#' blocks and main plots.
#'
#' @param x An `"ad_layout"` object from [make_design()].
#' @param color_by Column of the field book that colours the plots
#'   (default: `A` in split and factorial designs, otherwise `Treatment`).
#' @param labels Write the treatment in each plot.
#' @param numbers Write the plot number in each plot.
#' @param main Title.
#' @param ... Not used.
#' @return Invisibly, the field book.
#' @examples
#' plot(make_design("strip", A = c("Drip", "Furrow", "Sprinkler"), B = c("V1", "V2"), r = 2, seed = 3))
#' @export
plot.ad_layout <- function(x, color_by = NULL, labels = TRUE, numbers = TRUE, main = NULL, ...) {
  b <- x$book
  if (is.null(color_by)) color_by <- if ("A" %in% names(b) && x$design %in% c("split_rcbd", "strip", "splitsplit")) "A" else "Treatment"
  lv <- unique(b[[color_by]]); lv <- lv[order(match(lv, unique(b[[color_by]][order(b$Plot)])))]
  pal <- if (length(lv) <= 10) agri_colors(length(lv)) else grDevices::hcl.colors(length(lv), "Dark 3")
  fill <- .lighten(pal[match(b[[color_by]], lv)], 0.45)
  if ("Type" %in% names(b)) fill[b$Type == "entry"] <- "#f4f1ea"
  nr <- x$size[["rows"]]; nc <- x$size[["cols"]]
  op <- .style_par(mar = c(3.2, 3, 3, 1)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, nc + 0.5), ylim = c(nr + 0.5, 0.5), xaxs = "i", yaxs = "i", asp = if (nc / nr < 4) 1 else NA)
  graphics::rect(b$FieldCol - 0.47, b$FieldRow - 0.47, b$FieldCol + 0.47, b$FieldRow + 0.47, col = fill, border = "#ffffff")
  cx <- min(1, 9 / max(nchar(as.character(b$Treatment))))
  if (labels) graphics::text(b$FieldCol, b$FieldRow + if (numbers) 0.12 else 0,
                             if (x$design %in% c("split_rcbd", "splitsplit")) b[[if (x$design == "split_rcbd") "B" else "C"]]
                             else if (x$design == "graeco") paste(b$Treatment, b$Greek) else b$Treatment,
                             cex = 0.8 * cx, font = 2, col = "#1b1b1b")
  if (numbers) graphics::text(b$FieldCol - 0.4, b$FieldRow - 0.4, b$Plot, adj = c(0, 1), cex = 0.6, col = "#5f6368")
  bx <- x$boxes
  if (!is.null(bx)) {
    m <- bx[bx$kind != "block", ]
    if (nrow(m)) graphics::rect(m$c0 - 0.5, m$r0 - 0.5, m$c1 + 0.5, m$r1 + 0.5, border = "#3c4043",
                                lwd = ifelse(m$kind == "main", 1.6, 0.8), lty = ifelse(m$kind == "main", 1, 2))
    k <- bx[bx$kind == "block", ]
    if (nrow(k) > 1 || x$design != "crd") {
      graphics::rect(k$c0 - 0.5, k$r0 - 0.5, k$c1 + 0.5, k$r1 + 0.5, border = "#1b1b1b", lwd = 2.4)
      # stacked blocks touch each other: their names go in the left margin
      if (all(k$c0 == 1)) graphics::text(0.38, (k$r0 + k$r1) / 2, k$name, srt = 90, adj = c(0.5, 0), cex = 0.75, font = 2, xpd = NA)
      else graphics::text(k$c0 - 0.45, k$r0 - 0.5, k$name, adj = c(0, -0.4), cex = 0.75, font = 2, xpd = NA)
    }
  }
  graphics::title(main = if (is.null(main)) sprintf("%s, seed %s", sub(" \\(.*", "", x$design_name), x$seed) else main,
                  font.main = 1, cex.main = 1.05, line = 1.6)
  if (color_by != "Treatment" || length(lv) <= 12)
    graphics::legend("bottom", legend = lv, fill = .lighten(pal, 0.45), border = NA, horiz = TRUE, bty = "n",
                     cex = 0.75, xpd = NA, inset = c(0, -0.13), title = if (color_by == "Treatment") NULL else color_by)
  invisible(b)
}

# ---------- replicates and power ----------

.err_df <- function(design, t, r) if (design == "rcbd") (t - 1) * (r - 1) else t * (r - 1)

#' Replicates and power for a difference between two treatments
#'
#' Power of the two-sided t-test to detect a difference of `diff` per cent
#' of the mean between two treatment means, when the coefficient of
#' variation of the trial is `cv` per cent, as a function of the number of
#' replicates. The error degrees of freedom are those of a CRD or an RCBD
#' with `t` treatments, and the power uses the noncentral t distribution.
#'
#' @param cv Coefficient of variation expected (per cent), from earlier
#'   trials on the same site and crop.
#' @param diff Difference to detect, per cent of the mean.
#' @param t Number of treatments.
#' @param design `"rcbd"` or `"crd"`.
#' @param alpha Significance level.
#' @param power Power wanted.
#' @param r Replicates for which the table is computed.
#' @return An object of class `"ad_power"` with the replicates needed
#'   (`replicates`, `NA` above 200) and a table with the error df, the power
#'   and the detectable difference for each number of replicates.
#' @examples
#' power_design(cv = 15, diff = 15, t = 4)
#' plot(power_design(cv = 10, diff = 15, t = 4))
#' @export
power_design <- function(cv, diff, t, design = c("rcbd", "crd"), alpha = 0.05, power = 0.8, r = 2:30) {
  design <- match.arg(design)
  if (cv <= 0 || diff <= 0) stop("'cv' and 'diff' must be positive", call. = FALSE)
  pw <- function(rr) {
    df <- .err_df(design, t, rr); tc <- stats::qt(1 - alpha / 2, df); ncp <- diff / cv * sqrt(rr / 2)
    stats::pt(tc, df, ncp, lower.tail = FALSE) + stats::pt(-tc, df, ncp)
  }
  det <- function(rr) {
    df <- .err_df(design, t, rr)
    stats::uniroot(function(d) { tc <- stats::qt(1 - alpha / 2, df); ncp <- d / cv * sqrt(rr / 2)
      stats::pt(tc, df, ncp, lower.tail = FALSE) + stats::pt(-tc, df, ncp) - power }, c(1e-6, 100 * cv), tol = 1e-10)$root
  }
  need <- NA_integer_
  for (rr in 2:200) if (pw(rr) >= power) { need <- rr; break }
  tab <- data.frame(r = r, df_error = .err_df(design, t, r), power = vapply(r, pw, 0), detectable = vapply(r, det, 0))
  structure(list(replicates = need, table = tab, cv = cv, diff = diff, t = t, design = design, alpha = alpha, target = power),
            class = "ad_power")
}

#' @export
print.ad_power <- function(x, ...) {
  cat(sprintf("\nReplicates to detect a difference of %g %% of the mean (CV %g %%, %d treatments, %s, alpha %g, power %g)\n",
              x$diff, x$cv, x$t, toupper(x$design), x$alpha, x$target))
  cat(if (is.na(x$replicates)) "More than 200 replicates would be needed: the difference is too small for that CV.\n"
      else sprintf("Replicates needed: %d (error df %d)\n", x$replicates, .err_df(x$design, x$t, x$replicates)))
  tab <- x$table; tab$power <- round(tab$power, 3); tab$detectable <- round(tab$detectable, 1)
  names(tab) <- c("Replicates", "Error df", "Power", "Detectable diff. (%)")
  cat("\n"); print(utils::head(tab, 15), row.names = FALSE)
  cat("\nHalving the CV is worth as much as quadrupling the replicates.\n")
  invisible(x)
}

#' @export
plot.ad_power <- function(x, main = NULL, ...) {
  op <- .style_par(mfrow = c(1, 2)); on.exit(graphics::par(op))
  tab <- x$table
  graphics::plot(tab$r, tab$power, type = "n", ylim = c(0, 1), xlab = "Replicates", ylab = "Power",
                 main = if (is.null(main)) sprintf("Power for %g %% (CV %g %%)", x$diff, x$cv) else main, font.main = 1, ...)
  graphics::abline(h = x$target, lty = 2, col = "#5f6368")
  graphics::lines(tab$r, tab$power, lwd = 2.2, col = "#2f7d4f")
  graphics::points(tab$r, tab$power, pch = 21, bg = "#2f7d4f", col = "#ffffff")
  if (!is.na(x$replicates)) graphics::abline(v = x$replicates, lty = 3, col = "#b5432f")
  graphics::plot(tab$r, tab$detectable, type = "n", xlab = "Replicates", ylab = "Detectable difference (% of mean)",
                 main = sprintf("Detectable difference (power %g)", x$target), font.main = 1)
  graphics::abline(h = x$diff, lty = 2, col = "#5f6368")
  graphics::lines(tab$r, tab$detectable, lwd = 2.2, col = "#c8842a")
  graphics::points(tab$r, tab$detectable, pch = 21, bg = "#c8842a", col = "#ffffff")
  invisible(tab)
}
