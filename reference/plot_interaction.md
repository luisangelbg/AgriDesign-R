# Interaction plot of two treatment factors

Least-squares means of every combination of `x` (horizontal axis) and
`trace` (one line per level), with error bars.

## Usage

``` r
plot_interaction(
  fit,
  x,
  trace,
  error = c("se", "ci", "none"),
  col = NULL,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL
)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- x, trace:

  Treatment factors on the horizontal axis and as lines.

- error:

  Error bars: `"se"` (standard error of the mean), `"ci"` (95 %
  confidence interval with the error degrees of freedom) or `"none"`.

- col:

  Fill colours, one per level (default:
  [`agri_colors()`](https://luisangelbg.github.io/AgriDesign-R/reference/agri_colors.md)).

- main, xlab, ylab:

  Titles; `NULL` builds them from the analysis.

- ylim:

  Range of the vertical axis.

## Value

Invisibly, the cell means plotted.

## Examples

``` r
fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
                    block = "Block")
plot_interaction(fit, "Fertilizer", "Variety")
```
