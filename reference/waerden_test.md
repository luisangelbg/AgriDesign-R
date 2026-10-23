# van der Waerden normal-scores test

Compares the treatments of a one-way layout (no blocks) after replacing
every observation by its normal score `qnorm(R / (N + 1))`, where `R` is
its rank among all `N` observations (average ranks for ties). The
statistic is chi-square with `t - 1` degrees of freedom; the pairwise
comparisons are Conover's t tests on the mean scores with `N - t`
degrees of freedom. It is as powerful as the F test under normality and
more robust than it.

## Usage

``` r
waerden_test(
  x,
  response = NULL,
  treatment = NULL,
  alpha = 0.05,
  p_adjust = "holm"
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

## Value

An object of class `c("ad_waerden", "ad_nonparametric")`; the groups
carry the mean normal score of every level (`mean_score`).

## References

van der Waerden, B. L. (1952). Order tests for the two-sample problem
and their power. *Indagationes Mathematicae* 14, 453-458.

Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd ed.
Wiley, New York, section 5.10.

## Examples

``` r
waerden_test(bean_varieties, "Yield_kg_plot", "Variety")
#> 
#> van der Waerden normal-scores test
#> Response: Yield_kg_plot   factor: Variety
#> chisq = 17.07, df = 4, p = 0.001874
#> 
#> Pairwise comparisons: Conover t on mean normal scores (p adjusted by holm)
#>         level_i        level_j statistic         p     p_adj different
#>    Negro Jamapa   Flor de Mayo    -4.675 2.989e-04 1.495e-03      TRUE
#>    Negro Jamapa Pinto Saltillo    -2.542 2.257e-02 4.514e-02      TRUE
#>    Negro Jamapa           Bayo     2.855 1.206e-02 3.659e-02      TRUE
#>    Negro Jamapa       Azufrado    -7.666 1.451e-06 1.305e-05      TRUE
#>    Flor de Mayo Pinto Saltillo     2.134 4.978e-02 4.978e-02      TRUE
#>    Flor de Mayo           Bayo     7.530 1.801e-06 1.441e-05      TRUE
#>    Flor de Mayo       Azufrado    -2.990 9.148e-03 3.659e-02      TRUE
#>  Pinto Saltillo           Bayo     5.396 7.423e-05 5.196e-04      TRUE
#>  Pinto Saltillo       Azufrado    -5.124 1.246e-04 7.477e-04      TRUE
#>            Bayo       Azufrado   -10.520 2.553e-08 2.553e-07      TRUE
#> 
#> Groups (levels sharing a letter do not differ at alpha = 0.05):
#>           level n median  mean mean_rank mean_score group
#>        Azufrado 4  2.820 2.812     18.25    1.18939     a
#>    Flor de Mayo 4  2.595 2.598     14.50    0.51321     b
#>  Pinto Saltillo 4  2.425 2.420     10.75    0.03074     c
#>    Negro Jamapa 4  2.115 2.132      6.25   -0.54395     d
#>            Bayo 4  1.920 1.917      2.75   -1.18939     e
#> 
#> Variance of the normal scores = 0.7943; pairwise t with 15 df
```
