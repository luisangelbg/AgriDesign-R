# Plot the pairwise differences with their intervals

Every pairwise difference of a comparison with its interval: the
critical difference of the procedure (simultaneous for Tukey, Scheffe,
Dunnett...) or, for Holm and Games-Howell, the 95 % t interval.
Differences declared significant are drawn in red.

## Usage

``` r
plot_differences(
  x,
  sort = TRUE,
  main = NULL,
  xlab = "Difference between means",
  col = c("#2f7d4f", "#b5432f")
)
```

## Arguments

- x:

  An `"ad_comparison"` object from
  [`compare_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/compare_means.md).

- sort:

  Sort by the size of the difference.

- main, xlab:

  Titles.

- col:

  Colours of non-significant and significant differences.

## Value

Invisibly, the differences and interval limits plotted.

## Examples

``` r
fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
plot_differences(compare_means(fit, "Variety"))
```
