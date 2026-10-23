# Catalogue of experimental designs.
# Every design knows how many treatment factors and block columns it needs and
# builds its model terms from the roles. Each term carries the error stratum
# that tests it (`error`) and whether it is itself an error stratum (`is_error`).

.term <- function(factors, error = "Residuals", is_error = FALSE, label = NULL,
                  numeric = character(0)) {
  name <- paste(c(factors, numeric), collapse = ":")
  list(name = name, factors = factors, numeric = numeric, error = error,
       is_error = is_error, label = if (is.null(label)) name else label)
}

.cov_terms <- function(r) {
  lapply(r$covariate, function(x) .term(character(0), numeric = x))
}

.factorial_terms <- function(f) {
  out <- list()
  for (k in seq_along(f)) {
    cmb <- utils::combn(f, k, simplify = FALSE)
    out <- c(out, lapply(cmb, .term))
  }
  out
}

.catalog <- list(
  crd = list(
    name = "Completely randomised design (CRD)",
    factors = c(1, 1), blocks = c(0, 0),
    build = function(r, o) c(.cov_terms(r), list(.term(r$treatment[1])))),
  rcbd = list(
    name = "Randomised complete block design (RCBD)",
    factors = c(1, 1), blocks = c(1, 1),
    build = function(r, o) c(.cov_terms(r), list(.term(r$block[1]), .term(r$treatment[1])))),
  grcbd = list(
    name = "Generalised RCBD (several plots per block x treatment)",
    factors = c(1, 1), blocks = c(1, 1), need_reps = TRUE,
    build = function(r, o) {
      B <- r$block[1]; A <- r$treatment[1]
      BA <- paste(B, A, sep = ":")
      c(.cov_terms(r), list(
        .term(B),
        .term(A, error = if (o$random_blocks) BA else "Residuals"),
        .term(c(B, A), is_error = o$random_blocks,
              label = if (o$random_blocks) "Block x treatment (Error a)" else "Block x treatment")))
    }),
  latin = list(
    name = "Latin square",
    factors = c(1, 1), blocks = c(0, 0), need_row_col = TRUE,
    build = function(r, o) c(.cov_terms(r), list(.term(r$row), .term(r$col), .term(r$treatment[1])))),
  graeco = list(
    name = "Graeco-Latin square",
    factors = c(1, 1), blocks = c(1, 1), need_row_col = TRUE,
    build = function(r, o) c(.cov_terms(r), list(.term(r$row), .term(r$col), .term(r$block[1]), .term(r$treatment[1])))),
  youden = list(
    name = "Youden square (incomplete Latin square), intra-block analysis",
    factors = c(1, 1), blocks = c(0, 0), need_row_col = TRUE,
    build = function(r, o) c(.cov_terms(r), list(.term(r$row), .term(r$col), .term(r$treatment[1])))),
  bib = list(
    name = "Incomplete blocks without replicate grouping (BIBD, cyclic), intra-block analysis",
    factors = c(1, 1), blocks = c(1, 1),
    build = function(r, o) c(.cov_terms(r), list(.term(r$block[1]), .term(r$treatment[1])))),
  ibd = list(
    name = "Incomplete block design (lattice, alpha, BIBD), intra-block analysis",
    factors = c(1, 1), blocks = c(2, 2),
    build = function(r, o) {
      R <- r$block[1]; B <- r$block[2]
      c(.cov_terms(r), list(.term(R), .term(c(R, B), label = paste0(B, "(", R, ")")),
                            .term(r$treatment[1])))
    }),
  augmented = list(
    name = "Augmented design (replicated checks + unreplicated entries)",
    factors = c(1, 1), blocks = c(1, 1),
    build = function(r, o) c(.cov_terms(r), list(.term(r$block[1]), .term(r$treatment[1])))),
  fact_crd = list(
    name = "Factorial in CRD",
    factors = c(2, 4), blocks = c(0, 0),
    build = function(r, o) c(.cov_terms(r), .factorial_terms(r$treatment))),
  fact_rcbd = list(
    name = "Factorial in RCBD",
    factors = c(2, 4), blocks = c(1, 1),
    build = function(r, o) c(.cov_terms(r), list(.term(r$block[1])), .factorial_terms(r$treatment))),
  split_rcbd = list(
    name = "Split-plot in RCBD (main plots in blocks)",
    factors = c(2, 2), blocks = c(1, 1),
    build = function(r, o) {
      A <- r$treatment[1]; B <- r$treatment[2]; R <- r$block[1]
      ea <- paste(R, A, sep = ":")
      c(.cov_terms(r), list(
        .term(R, error = ea), .term(A, error = ea),
        .term(c(R, A), is_error = TRUE, label = "Error a (main plot)"),
        .term(B), .term(c(A, B))))
    }),
  split_crd = list(
    name = "Split-plot in CRD (main plots completely randomised)",
    factors = c(2, 2), blocks = c(1, 1),
    build = function(r, o) {
      A <- r$treatment[1]; B <- r$treatment[2]; R <- r$block[1]
      ea <- paste(A, R, sep = ":")
      c(.cov_terms(r), list(
        .term(A, error = ea),
        .term(c(A, R), is_error = TRUE, label = paste0("Error a (replicates within ", A, ")")),
        .term(B), .term(c(A, B))))
    }),
  strip_rcbd = list(
    name = "Strip-plot (split-block) in RCBD",
    factors = c(2, 2), blocks = c(1, 1),
    build = function(r, o) {
      A <- r$treatment[1]; B <- r$treatment[2]; R <- r$block[1]
      ea <- paste(R, A, sep = ":"); eb <- paste(R, B, sep = ":")
      c(.cov_terms(r), list(
        .term(R), .term(A, error = ea), .term(c(R, A), is_error = TRUE, label = "Error a"),
        .term(B, error = eb), .term(c(R, B), is_error = TRUE, label = "Error b"),
        .term(c(A, B))))
    }),
  splitsplit_rcbd = list(
    name = "Split-split-plot in RCBD",
    factors = c(3, 3), blocks = c(1, 1),
    build = function(r, o) {
      A <- r$treatment[1]; B <- r$treatment[2]; C <- r$treatment[3]; R <- r$block[1]
      ea <- paste(R, A, sep = ":"); eb <- paste(R, A, B, sep = ":")
      c(.cov_terms(r), list(
        .term(R, error = ea), .term(A, error = ea),
        .term(c(R, A), is_error = TRUE, label = "Error a"),
        .term(B, error = eb), .term(c(A, B), error = eb),
        .term(c(R, A, B), is_error = TRUE, label = "Error b"),
        .term(C), .term(c(A, C)), .term(c(B, C)), .term(c(A, B, C))))
    }),
  nested = list(
    name = "Nested (hierarchical): second factor within the first",
    factors = c(2, 2), blocks = c(0, 0),
    build = function(r, o) {
      A <- r$treatment[1]; B <- r$treatment[2]
      c(.cov_terms(r), list(
        .term(A, error = if (o$random_nested) paste(A, B, sep = ":") else "Residuals"),
        .term(c(A, B), is_error = o$random_nested, label = paste0(B, "(", A, ")"))))
    })
)

#' List the designs that can be analysed or randomised
#'
#' `designs()` lists the 16 designs that [design_anova()] analyses;
#' `designs("layout")` lists the 16 designs that [make_design()] randomises,
#' with the arguments each one needs and the design to give to
#' [design_anova()] once the responses are in the field book.
#'
#' Designs analysed by [design_anova()] (`design =`):
#'
#' | id | Design | Roles |
#' |----|--------|-------|
#' | `crd` | Completely randomised | treatment |
#' | `rcbd` | Randomised complete blocks | treatment, block |
#' | `grcbd` | Generalised RCBD (replicates within blocks) | treatment, block |
#' | `latin` | Latin square | treatment, row, col |
#' | `graeco` | Graeco-Latin square | treatment, row, col, block (Greek letters) |
#' | `youden` | Youden square | treatment, row, col |
#' | `bib` | Incomplete blocks (BIBD, cyclic) | treatment, block |
#' | `ibd` | Resolvable incomplete blocks (lattice, alpha) | treatment, block = c(rep, block) |
#' | `augmented` | Augmented (checks and new entries) | treatment, block |
#' | `fact_crd` | Factorial in CRD (2 to 4 factors) | treatment = c(A, B, ...) |
#' | `fact_rcbd` | Factorial in RCBD (2 to 4 factors) | treatment = c(A, B, ...), block |
#' | `split_rcbd` | Split-plot in RCBD | treatment = c(main, sub), block |
#' | `split_crd` | Split-plot in CRD | treatment = c(main, sub), block (replicate) |
#' | `strip_rcbd` | Strip-plot (split-block) | treatment = c(A, B), block |
#' | `splitsplit_rcbd` | Split-split-plot in RCBD | treatment = c(main, sub, sub-sub), block |
#' | `nested` | Nested (hierarchical) | treatment = c(outer, inner) |
#'
#' Designs randomised by [make_design()] (first argument; `C` is optional in
#' factorials, and the conditions on the number of treatments and the block
#' size are in the `conditions` column of `designs("layout")`):
#'
#' | id | Design | Arguments | Analyse as |
#' |----|--------|-----------|------------|
#' | `crd` | Completely randomised | `trt`, `r` | `crd` |
#' | `rcbd` | Randomised complete blocks | `trt`, `r` | `rcbd` |
#' | `latin` | Latin square | `trt` | `latin` |
#' | `graeco` | Graeco-Latin square | `trt`, `greek` | `graeco` |
#' | `youden` | Youden square | `trt`, `k` | `youden` |
#' | `bib` | Balanced incomplete blocks | `trt`, `k` | `bib` |
#' | `cyclic` | Cyclic incomplete blocks | `trt`, `k`, `r` | `bib` |
#' | `ibd` | Resolvable incomplete blocks | `trt`, `k`, `r` | `ibd` |
#' | `alpha` | Alpha design | `trt`, `k`, `r` | `ibd` |
#' | `lattice` | Square lattice | `trt`, `k`, `r` | `ibd` |
#' | `augmented` | Augmented | `checks`, `entries`, `r` | `augmented` |
#' | `fact_crd` | Factorial in CRD | `A`, `B`, `C`, `r` | `fact_crd` |
#' | `fact_rcbd` | Factorial in RCBD | `A`, `B`, `C`, `r` | `fact_rcbd` |
#' | `split_rcbd` | Split-plot in RCBD | `A`, `B`, `r` | `split_rcbd` |
#' | `strip` | Strip-plot | `A`, `B`, `r` | `strip_rcbd` |
#' | `splitsplit` | Split-split-plot | `A`, `B`, `C`, `r` | `splitsplit_rcbd` |
#'
#' @param type `"analysis"` (designs of [design_anova()]) or `"layout"`
#'   (designs of [make_design()]).
#' @return A data frame. For `"analysis"`: the identifier of each design (the
#'   value to pass to `design =` in [design_anova()]), its name, the number of
#'   treatment factors and block columns, and whether it needs rows and
#'   columns. For `"layout"`: the identifier (the first argument of
#'   [make_design()]), its name, the arguments it needs (in parentheses when
#'   optional), its conditions and the design of [design_anova()] that
#'   analyses its field book.
#' @examples
#' designs()
#' designs("layout")
#'
#' # the layout reports the call that analyses it
#' make_design("strip", A = c("Drip", "Furrow"), B = c("V1", "V2", "V3"), r = 3, seed = 1)$analysis
#' @export
designs <- function(type = c("analysis", "layout")) {
  if (match.arg(type) == "layout") {
    id <- names(.generated)
    return(data.frame(id = id, name = unname(.generated), needs = unname(.layout_needs[id]),
                      conditions = unname(.layout_conditions[id]),
                      analyse_as = vapply(id, function(d) sub('.*design = "([a-z_]+)".*', "\\1", .analysis_call(d)), ""),
                      row.names = NULL, stringsAsFactors = FALSE))
  }
  data.frame(
    id = names(.catalog),
    name = vapply(.catalog, `[[`, "", "name"),
    treatments = vapply(.catalog, function(d) {
      if (d$factors[1] == d$factors[2]) as.character(d$factors[1])
      else paste(d$factors, collapse = "-")
    }, ""),
    blocks = vapply(.catalog, function(d) d$blocks[1], 0),
    row_col = vapply(.catalog, function(d) isTRUE(d$need_row_col), NA),
    row.names = NULL, stringsAsFactors = FALSE)
}

# Which designs fit the roles: returns a data frame (id, ok, why) and the suggestion.
.applicable <- function(r, data) {
  nf <- length(r$treatment); nb <- length(r$block)
  multi <- FALSE
  if (nf && nb) {
    key <- interaction(data[c(r$block[1], r$treatment)], drop = TRUE)
    multi <- any(table(key) > 1)
  }
  chk <- lapply(names(.catalog), function(id) {
    d <- .catalog[[id]]; why <- character(0)
    if (nf < d$factors[1] || nf > d$factors[2])
      why <- c(why, sprintf("needs %s treatment factor(s)",
                            if (d$factors[1] == d$factors[2]) d$factors[1] else paste(d$factors, collapse = "-")))
    if (nb < d$blocks[1]) why <- c(why, sprintf("needs %d block column(s)", d$blocks[1]))
    if (isTRUE(d$need_row_col) && (is.null(r$row) || is.null(r$col))) why <- c(why, "needs row and column")
    if (isTRUE(d$need_reps) && nb >= d$blocks[1] && !multi) why <- c(why, "needs several plots per block x treatment")
    data.frame(id = id, ok = !length(why), why = paste(why, collapse = "; "), stringsAsFactors = FALSE)
  })
  incomplete <- function(b) nf == 1 && any(table(data[[b]], data[[r$treatment[1]]]) == 0)
  sug <- if (!is.null(r$row) && !is.null(r$col) && nf == 1) {
           if (nb == 1) "graeco" else if (incomplete(r$col)) "youden" else "latin"
         }
         else if (nf == 1) { if (nb == 0) "crd" else if (nb >= 2) "ibd" else if (multi) "grcbd"
                             else if (incomplete(r$block[1])) "bib" else "rcbd" }
         else if (nb) "fact_rcbd" else "fact_crd"
  list(table = do.call(rbind, chk), suggested = sug)
}
