# Further designs for make_design(): Graeco-Latin squares, Youden squares,
# balanced incomplete blocks, cyclic designs, alpha designs and square
# lattices; and the efficiency factor of any block design.

# ---------- finite fields GF(2^m), for Graeco-Latin squares of order 4 and 8 ----------
.gf2_mult <- function(a, b, m, poly) {
  p <- 0L
  for (i in seq_len(m)) {
    if (bitwAnd(b, 1L)) p <- bitwXor(p, a)
    b <- bitwShiftR(b, 1L)
    a <- bitwShiftL(a, 1L)
    if (bitwAnd(a, bitwShiftL(1L, m))) a <- bitwXor(a, poly)
  }
  p
}

# Two orthogonal Latin squares of order t (symbols 1..t): over Z_t with
# L1 = i + j, L2 = i + 2j for odd t, and over GF(4) or GF(8) with
# L1 = i + j, L2 = i + w j (w a primitive element) for t = 4 and 8.
.mols_pair <- function(t) {
  i <- rep(0:(t - 1), each = t); j <- rep(0:(t - 1), t)
  if (t %% 2 == 1) {
    L1 <- (i + j) %% t; L2 <- (i + 2 * j) %% t
  } else if (t %in% c(4, 8)) {
    m <- log2(t); poly <- if (t == 4) 7L else 11L          # x^2 + x + 1, x^3 + x + 1
    L1 <- bitwXor(i, j)
    L2 <- bitwXor(i, vapply(j, function(v) .gf2_mult(2L, v, m, poly), 0L))
  } else stop("Graeco-Latin squares are built here for odd orders and for 4 and 8 ",
              "(none exists of order 2 or 6)", call. = FALSE)
  list(L1 = matrix(L1 + 1, t, t, byrow = TRUE), L2 = matrix(L2 + 1, t, t, byrow = TRUE))
}

# ---------- cyclic difference sets (v, k, lambda) for Youden squares and BIBDs ----------
.difference_sets <- list(
  "7-3" = c(0, 1, 3), "7-4" = c(0, 1, 2, 4), "11-5" = c(1, 3, 4, 5, 9), "11-6" = c(0, 2, 6, 7, 8, 10),
  "13-4" = c(0, 1, 3, 9), "13-9" = c(2, 4, 5, 6, 7, 8, 10, 11, 12), "15-7" = c(0, 1, 2, 4, 5, 8, 10),
  "15-8" = c(3, 6, 7, 9, 11, 12, 13, 14), "19-9" = c(1, 4, 5, 6, 7, 9, 11, 16, 17),
  "21-5" = c(3, 6, 7, 12, 14), "23-11" = c(1, 2, 3, 4, 6, 8, 9, 12, 13, 16, 18),
  "31-6" = c(1, 5, 11, 24, 25, 27))

.find_difference_set <- function(t, k) {
  if (k == t - 1) return(seq_len(t - 1) - 1)
  .difference_sets[[paste(t, k, sep = "-")]]
}

# ---------- efficiency factor ----------
# Canonical efficiency factors: non-zero eigenvalues of R^(-1/2) C R^(-1/2), with
# C = R - N K^-1 N' the information matrix of the treatments adjusted for blocks.
.efficiency <- function(trt, block) {
  trt <- factor(trt); block <- factor(block)
  N <- unclass(table(trt, block))
  rr <- rowSums(N); kk <- colSums(N)
  C <- diag(rr) - N %*% diag(1 / kk, length(kk)) %*% t(N)
  Rh <- diag(1 / sqrt(rr))
  e <- sort(eigen(Rh %*% C %*% Rh, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  e <- e[seq_len(nlevels(trt) - 1)]
  if (any(e < 1e-9)) return(list(E = 0, factors = e, connected = FALSE))
  list(E = (nlevels(trt) - 1) / sum(1 / e), factors = e, connected = TRUE)
}

#' Efficiency factor of a block design
#'
#' The efficiency factor of a block design is the harmonic mean of its
#' canonical efficiency factors: the average variance of a treatment
#' difference in a complete block design with the same replication, divided
#' by the average variance in this design (intra-block analysis, same error
#' variance). It is 1 for complete blocks and the target of incomplete block
#' designs. For resolvable designs the upper bound of Patterson and Williams
#' (1976) is also reported.
#'
#' @param x An `"ad_layout"` object from [make_design()], or a data frame.
#' @param treatment,block Column names of the treatment and of the
#'   (incomplete) block when `x` is a data frame. For a layout the blocks are
#'   the incomplete blocks when there are any.
#' @return A list with the efficiency factor `E`, the canonical efficiency
#'   factors, whether the design is connected, and `upper_bound` for
#'   equireplicated resolvable designs.
#' @examples
#' lay <- make_design("alpha", trt = 12, r = 2, k = 3, seed = 1)
#' design_efficiency(lay)$E
#' @export
design_efficiency <- function(x, treatment = "Treatment", block = NULL) {
  if (inherits(x, "ad_layout")) {
    b <- x$book
    blk <- if ("IncBlock" %in% names(b)) paste(b$Block, b$IncBlock) else if ("Block" %in% names(b)) b$Block
           else if ("Column" %in% names(b)) b$Column else stop("this layout has no blocks", call. = FALSE)
    trt <- b$Treatment
    reps <- if ("IncBlock" %in% names(b)) length(unique(b$Block)) else NA
  } else {
    .check_columns(x, c(treatment, block))
    trt <- x[[treatment]]; blk <- x[[block]]; reps <- NA
  }
  e <- .efficiency(trt, blk)
  t <- length(unique(trt)); k <- stats::median(table(blk))
  if (!is.na(reps)) {
    s <- t / k
    e$upper_bound <- (t - 1) * (reps - 1) / ((t - 1) * (reps - 1) + reps * (s - 1))
  }
  e
}

# ---------- generators (same format as .gen in generator.R) ----------
.gen_extra <- list(
  graeco = function(o) {
    t <- length(o$trt)
    if (length(o$greek) != t) stop("give as many Greek letters ('greek') as treatments", call. = FALSE)
    P <- .mols_pair(t)
    rp <- .shuffle(seq_len(t)); cp <- .shuffle(seq_len(t)); sl <- .shuffle(o$trt); sg <- .shuffle(o$greek)
    g <- expand.grid(j = seq_len(t), i = seq_len(t))
    cells <- .cells(sl[P$L1[cbind(rp[g$i], cp[g$j])]], Greek = sg[P$L2[cbind(rp[g$i], cp[g$j])]],
                    Row = g$i, Column = g$j, rep = NA, block_index = 0, in_block = seq_len(nrow(g)) - 1)
    cells$r <- g$i - 1; cells$c <- g$j - 1
    list(list(name = "Square", cells = cells))
  },
  youden = function(o) {
    t <- length(o$trt); k <- o$k
    D <- .find_difference_set(t, k)
    if (is.null(D)) stop("no cyclic difference set is stored for ", t, " treatments in blocks of ", k,
                         "; available (t, k): ", paste(c(names(.difference_sets), "any (t, t-1)"), collapse = ", "), call. = FALSE)
    rp <- .shuffle(seq_len(k)); cp <- .shuffle(seq_len(t) - 1); sl <- .shuffle(o$trt)
    g <- expand.grid(j = seq_len(t), i = seq_len(k))
    sym <- (D[rp[g$i]] + cp[g$j]) %% t + 1
    cells <- .cells(sl[sym], Row = g$i, Column = g$j, rep = NA, block_index = 0, in_block = seq_len(nrow(g)) - 1)
    cells$r <- g$i - 1; cells$c <- g$j - 1
    list(list(name = "Youden square", cells = cells))
  },
  bib = function(o) {
    t <- length(o$trt); k <- o$k
    D <- .find_difference_set(t, k)
    blocks <- if (!is.null(D)) lapply(0:(t - 1), function(j) (D + j) %% t + 1)
              else if (choose(t, k) <= 500) utils::combn(t, k, simplify = FALSE)
              else stop("no balanced incomplete block design is built for ", t, " treatments in blocks of ", k,
                        " (", choose(t, k), " blocks would be needed); use design = \"alpha\" or \"cyclic\"", call. = FALSE)
    sl <- .shuffle(o$trt); blocks <- blocks[.shuffle(seq_along(blocks))]
    lapply(seq_along(blocks), function(b) {
      p <- sl[.shuffle(blocks[[b]])]
      list(name = paste("Block", b), cells = .to_grid(.cells(p, rep = b, block_index = b - 1, in_block = seq_along(p) - 1), length(p)))
    })
  },
  cyclic = function(o) {
    t <- length(o$trt); k <- o$k
    if (o$r %% k) stop("in a cyclic design 'r' must be a multiple of the block size k (each initial block gives r = k)", call. = FALSE)
    nib <- o$r / k
    # initial blocks contain 0; choose them, one by one, to maximise the efficiency factor
    cand <- if (choose(t - 1, k - 1) <= 3000) utils::combn(seq_len(t - 1), k - 1, simplify = FALSE)
            else lapply(seq_len(3000), function(i) sort(sample.int(t - 1, k - 1)))
    chosen <- list()
    develop <- function(B) unlist(lapply(B, function(b) lapply(0:(t - 1), function(j) (c(0, b) + j) %% t + 1)), recursive = FALSE)
    for (s in seq_len(nib)) {
      best <- -1; pick <- NULL
      for (cb in cand) {
        bl <- develop(c(chosen, list(cb)))
        e <- .efficiency(unlist(bl), rep(seq_along(bl), each = k))$E
        if (e > best + 1e-12) { best <- e; pick <- cb }
      }
      chosen <- c(chosen, list(pick))
    }
    blocks <- develop(chosen)
    sl <- .shuffle(o$trt); blocks <- blocks[.shuffle(seq_along(blocks))]
    lapply(seq_along(blocks), function(b) {
      p <- sl[.shuffle(blocks[[b]])]
      list(name = paste("Block", b), cells = .to_grid(.cells(p, rep = b, block_index = b - 1, in_block = seq_along(p) - 1), length(p)))
    })
  },
  alpha = function(o) {
    t <- length(o$trt); k <- o$k; r <- o$r
    if (t %% k) stop("the number of treatments (", t, ") must be a multiple of the block size k (", k, ")", call. = FALSE)
    s <- t / k
    if (k > s) stop("alpha designs need k <= t / k (blocks no larger than their number in a replicate)", call. = FALSE)
    blocks_of <- function(A) unlist(lapply(seq_len(r), function(j) lapply(0:(s - 1), function(l)
      (seq_len(k) - 1) * s + (A[, j] + l) %% s + 1)), recursive = FALSE)
    eff <- function(A) { bl <- blocks_of(A); .efficiency(unlist(bl), rep(seq_along(bl), each = k))$E }
    # reduced alpha array: first row and first column zero; the rest searched for the highest efficiency
    if (isTRUE(o$lattice)) {
      A <- outer(0:(k - 1), 0:(r - 1)) %% s
    } else {
      m <- (k - 1) * (r - 1)
      to_A <- function(v) { Ac <- matrix(0, k, r); Ac[-1, -1] <- v; Ac }
      if (s^m <= 2000) {
        # small problems: every reduced alpha array
        grid <- as.matrix(expand.grid(rep(list(0:(s - 1)), m)))
        es <- apply(grid, 1, function(v) eff(to_A(v)))
        A <- to_A(grid[which.max(es), ])
      } else {
        # random arrays, then coordinate ascent (one entry at a time) from the five best
        starts <- lapply(seq_len(300), function(i) sample.int(s, m, replace = TRUE) - 1)
        es <- vapply(starts, function(v) eff(to_A(v)), 0)
        best <- -1
        for (v in starts[order(-es)[1:5]]) {
          cur <- eff(to_A(v))
          repeat {
            improved <- FALSE
            for (pos in seq_len(m)) for (val in 0:(s - 1)) {
              if (val == v[pos]) next
              w <- v; w[pos] <- val; e <- eff(to_A(w))
              if (e > cur + 1e-12) { v <- w; cur <- e; improved <- TRUE }
            }
            if (!improved) break
          }
          if (cur > best + 1e-12) { best <- cur; A <- to_A(v) }
        }
      }
    }
    sl <- .shuffle(o$trt)
    lapply(seq_len(r), function(j) {
      bl <- lapply(0:(s - 1), function(l) (seq_len(k) - 1) * s + (A[, j] + l) %% s + 1)
      bl <- bl[.shuffle(seq_along(bl))]
      cells <- do.call(rbind, lapply(seq_along(bl), function(b) {
        p <- sl[.shuffle(bl[[b]])]
        .cells(p, rep = j, IncBlock = b, block_index = j - 1, r = b - 1, c = seq_len(k) - 1, in_block = (b - 1) * k + seq_len(k) - 1)
      }))
      list(name = paste("Rep", j), cells = cells,
           boxes = data.frame(name = paste0("B", seq_along(bl)), r0 = seq_along(bl) - 1, c0 = 0, r1 = seq_along(bl) - 1,
                              c1 = k - 1, kind = "sub", stringsAsFactors = FALSE))
    })
  })
