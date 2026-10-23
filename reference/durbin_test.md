# Durbin test for balanced incomplete blocks

Rank-based comparison of treatments in a balanced incomplete block
design (every block has `k` of the `t` treatments, every treatment
appears `r` times): the observations are ranked within each block and
the rank sums are compared (Durbin 1951). With complete blocks it is
Friedman's test. Pairwise comparisons use Conover's t test on the rank
sums with `b k - b - t + 1` degrees of freedom; ties get average ranks.

## Usage

``` r
durbin_test(
  x,
  response = NULL,
  treatment = NULL,
  block = NULL,
  alpha = 0.05,
  p_adjust = "holm"
)
```

## Arguments

- x:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
  with one treatment factor and blocks, or a data frame.

- response, treatment, block:

  Column names when `x` is a data frame; several block columns
  (replicate, block within replicate) are combined.

- alpha:

  Significance level of the pairwise comparisons.

- p_adjust:

  Adjustment of the pairwise p-values, as in
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); `"none"`
  gives Conover's protected comparisons (to be read only when the
  overall test is significant).

## Value

An object of class `c("ad_durbin", "ad_nonparametric")` with the test
(chi-square and its F form), the pairwise comparisons, the groups with
their rank sums and letters, the design parameters (`design`) and the
least significant difference of rank sums (`lsd_ranks`).

## References

Durbin, J. (1951). Incomplete blocks in ranking experiments. *British
Journal of Statistical Psychology* 4, 85-90.

Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd ed.
Wiley, New York, section 5.9.

## Examples

``` r
durbin_test(bib_sorghum, "Yield_t_ha", "Hybrid", "Block")
#> 
#> Durbin rank test for balanced incomplete blocks
#> Response: Yield_t_ha   factor: Hybrid
#> chisq = 10.44, df = 6, p = 0.1071
#> 
#> Pairwise comparisons: Conover t on rank sums (p adjusted by holm)
#>  level_i level_j statistic        p  p_adj different
#>       H1      H2   -0.3118 0.763160 1.0000     FALSE
#>       H1      H4   -1.8708 0.098282 1.0000     FALSE
#>       H1      H3    1.2472 0.247591 1.0000     FALSE
#>       H1      H5   -0.3118 0.763160 1.0000     FALSE
#>       H1      H6    1.8708 0.098282 1.0000     FALSE
#>       H1      H7   -0.6236 0.550247 1.0000     FALSE
#>       H2      H4   -1.5590 0.157609 1.0000     FALSE
#>       H2      H3    1.5590 0.157609 1.0000     FALSE
#>       H2      H5    0.0000 1.000000 1.0000     FALSE
#>       H2      H6    2.1826 0.060613 1.0000     FALSE
#>       H2      H7   -0.3118 0.763160 1.0000     FALSE
#>       H4      H3    3.1180 0.014271 0.2854     FALSE
#>       H4      H5    1.5590 0.157609 1.0000     FALSE
#>       H4      H6    3.7417 0.005692 0.1195     FALSE
#>       H4      H7    1.2472 0.247591 1.0000     FALSE
#>       H3      H5   -1.5590 0.157609 1.0000     FALSE
#>       H3      H6    0.6236 0.550247 1.0000     FALSE
#>       H3      H7   -1.8708 0.098282 1.0000     FALSE
#>       H5      H6    2.1826 0.060613 1.0000     FALSE
#>       H5      H7   -0.3118 0.763160 1.0000     FALSE
#>       H6      H7   -2.4944 0.037264 0.7080     FALSE
#> 
#> Groups (levels sharing a letter do not differ at alpha = 0.05):
#>  level n median  mean mean_rank rank_sum group
#>     H4 3   7.97 7.763     3.000      9.0     a
#>     H7 3   7.41 7.527     2.333      7.0     a
#>     H2 3   6.45 6.507     2.167      6.5     a
#>     H5 3   5.68 6.417     2.167      6.5     a
#>     H1 3   6.78 6.497     2.000      6.0     a
#>     H3 3   6.43 6.127     1.333      4.0     a
#>     H6 3   5.98 5.707     1.000      3.0     a
#> 
#> Design: 7 treatments in 7 blocks of 3 plots, r = 3, lambda = 1
#> F form: F = 3.917 on 6 and 8 df, p = 0.03973
```
