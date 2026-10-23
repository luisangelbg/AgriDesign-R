# Mood's median test

Compares the treatments of a one-way layout by the number of
observations above the grand median: a chi-square test on the 2 x t
table of counts above and at-or-below the median (Mood 1950; Conover
1999). Pairs of treatments are compared with the same test on their 2 x
2 table, either around the median of the two treatments together
(`pair_median = "pair"`, the default) or around the grand median
(`"grand"`), or with Fisher's exact test (`exact = TRUE`). It assumes
little more than independence, and resists outliers and unequal spreads,
at the cost of power.

## Usage

``` r
median_test(
  x,
  response = NULL,
  treatment = NULL,
  alpha = 0.05,
  p_adjust = "holm",
  correct = FALSE,
  pair_median = c("pair", "grand"),
  exact = FALSE
)
```

## Arguments

- x:

  An `"ad_anova"` object of a one-factor design without blocks, or a
  data frame.

- response, treatment:

  Column names when `x` is a data frame.

- alpha:

  Significance level of the pairwise comparisons.

- p_adjust:

  Adjustment of the pairwise p-values, as in
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); `"none"`
  gives Conover's protected comparisons (to be read only when the
  overall test is significant).

- correct:

  Yates' continuity correction in the 2 x 2 tables (the overall test is
  corrected only when there are two treatments).

- pair_median:

  Median used for the pairwise tables.

- exact:

  Use Fisher's exact test for the pairwise tables.

## Value

An object of class `c("ad_median_test", "ad_nonparametric")`; `table`
holds the counts above and at-or-below the grand median and
`grand_median` the median of all observations.

## References

Mood, A. M. (1950). *Introduction to the Theory of Statistics*.
McGraw-Hill, New York, 394-399.

Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd ed.
Wiley, New York, section 4.3.

## Examples

``` r
median_test(bean_varieties, "Yield_kg_plot", "Variety")
#> 
#> Mood's median test
#> Response: Yield_kg_plot   factor: Variety
#> chisq = 16, df = 4, p = 0.003019
#> 
#>           level n above at_or_below
#>    Negro Jamapa 4     0           4
#>    Flor de Mayo 4     4           0
#>  Pinto Saltillo 4     2           2
#>            Bayo 4     0           4
#>        Azufrado 4     4           0
#> 
#> Pairwise comparisons: median test by pairs (median of each pair) (p adjusted by holm)
#>         level_i        level_j statistic        p   p_adj different
#>    Negro Jamapa   Flor de Mayo         8 0.004678 0.04678      TRUE
#>    Negro Jamapa Pinto Saltillo         8 0.004678 0.04678      TRUE
#>    Negro Jamapa           Bayo         2 0.157299 0.47190     FALSE
#>    Negro Jamapa       Azufrado         8 0.004678 0.04678      TRUE
#>    Flor de Mayo Pinto Saltillo         2 0.157299 0.47190     FALSE
#>    Flor de Mayo           Bayo         8 0.004678 0.04678      TRUE
#>    Flor de Mayo       Azufrado         2 0.157299 0.47190     FALSE
#>  Pinto Saltillo           Bayo         8 0.004678 0.04678      TRUE
#>  Pinto Saltillo       Azufrado         8 0.004678 0.04678      TRUE
#>            Bayo       Azufrado         8 0.004678 0.04678      TRUE
#> 
#> Groups (levels sharing a letter do not differ at alpha = 0.05):
#>           level n median  mean mean_rank above group
#>        Azufrado 4  2.820 2.812     18.25     4     a
#>    Flor de Mayo 4  2.595 2.598     14.50     4    ab
#>  Pinto Saltillo 4  2.425 2.420     10.75     2     b
#>    Negro Jamapa 4  2.115 2.132      6.25     0     c
#>            Bayo 4  1.920 1.917      2.75     0     c
#> 
#> Grand median = 2.425; counts in 'table' are above / at-or-below it.
#> 
#> Some levels have fewer than 5 observations: chi-square p-values are rough; consider exact = TRUE.
```
