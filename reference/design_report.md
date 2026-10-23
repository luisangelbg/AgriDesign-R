# Write an HTML report of an analysis

Builds a self-contained HTML report of a fitted design: a materials and
methods paragraph ready to adapt, descriptive statistics, the assumption
checks with their grade, the analysis of variance with a draft of the
results, the means with their letters, trend analysis of quantitative
factors, simple effects of two-factor interactions, figures embedded as
vector graphics, appendices with the pairwise comparisons and the data,
and how to cite the software. Open it in any web browser and print it to
PDF from there (A4 page style included).

## Usage

``` r
design_report(
  fit,
  file = NULL,
  title = NULL,
  authors = "",
  affiliation = "",
  trial = "",
  objective = "",
  units = "",
  method = "tukey",
  alpha = 0.05,
  control = NULL,
  assumptions = TRUE,
  figures = TRUE,
  pairwise = TRUE,
  data = FALSE,
  cite = TRUE,
  open = FALSE
)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- file:

  Output file (`.html`); `NULL` only returns the HTML.

- title, authors, affiliation, trial, objective:

  Texts for the head of the report.

- units:

  Units of the response, for the methods paragraph.

- method, alpha:

  Mean-separation procedure and significance level (see
  [`compare_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/compare_means.md)).

- control:

  Control level for `method = "dunnett"`.

- assumptions, figures, pairwise, data, cite:

  Include those sections.

- open:

  Open the report in the browser when it is written.

## Value

The HTML text, invisibly.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
f <- tempfile(fileext = ".html")
design_report(fit, f, title = "Nitrogen rates in maize", units = "t/ha", figures = FALSE)
file.exists(f)
#> [1] TRUE
```
