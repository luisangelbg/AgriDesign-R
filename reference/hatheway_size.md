# Plot size needed to detect a difference (Hatheway)

Hatheway's (1961) relation between plot size, replication and the
difference that can be detected, from Smith's index:
`x^b = 2 (t1 + t2)^2 CV1^2 / (r d^2)`, where `CV1` is the coefficient of
variation of the basic units (%), `d` the true difference to detect as a
percentage of the mean, `r` the number of replicates, `t1` the two-sided
t for `alpha` and `t2` the one-sided t for the power wanted.

## Usage

``` r
hatheway_size(
  index,
  reps,
  difference,
  alpha = 0.05,
  power = 0.8,
  treatments = NULL,
  cv = NULL
)
```

## Arguments

- index:

  An `"ad_smith"` object from
  [`smith_index()`](https://luisangelbg.github.io/AgriDesign-R/reference/smith_index.md),
  or the value of `b` (then give `cv`).

- reps:

  Number(s) of replicates.

- difference:

  Difference(s) to detect, as a percentage of the mean.

- alpha:

  Significance level (two-sided).

- power:

  Probability of detecting the difference.

- treatments:

  Optional number of treatments: the t values then use the error degrees
  of freedom of a randomised complete block design,
  `(reps - 1) (treatments - 1)`; otherwise normal quantiles are used.

- cv:

  Coefficient of variation (%) of the basic units; by default the one of
  the uniformity trial.

## Value

A data frame with `reps`, `difference`, the error `df` and the plot size
in basic units (`plot_size`).

## References

Hatheway, W. H. (1961). Convenient plot size. *Agronomy Journal* 53,
279-280.

## Examples

``` r
sm <- smith_index(uniformity_trial)
hatheway_size(sm, reps = 3:6, difference = c(10, 15), treatments = 8)
#>   reps difference df  plot_size
#> 1    3         10 14 35.8278655
#> 2    4         10 21 15.4204465
#> 3    5         10 28  8.3095551
#> 4    6         10 35  5.0828259
#> 5    3         15 14  4.7268758
#> 6    4         15 21  2.0344649
#> 7    5         15 28  1.0963041
#> 8    6         15 35  0.6705922
```
