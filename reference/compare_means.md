# Compare treatment means

Separates the least-squares means of a treatment factor (or of the
combinations of two factors) with one of eleven procedures and labels
them with letters: levels that share a letter do not differ. Each factor
is compared with its own error stratum (for instance, main plots of a
split-plot with Error a).

## Usage

``` r
compare_means(
  fit,
  factor,
  method = "tukey",
  alpha = 0.05,
  control = NULL,
  by = NULL
)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- factor:

  Treatment factor whose means are compared; two names compare all their
  combinations together (not available for split- or strip-plots, where
  `by` should be used).

- method:

  Procedure, see Details.

- alpha:

  Significance level.

- control:

  Level used as control by Dunnett's test (default: first level).

- by:

  Optional second treatment factor: compare `factor` within each of its
  levels.

## Value

An object of class `"ad_comparison"` with the means and their letters
(`means`), the pairwise comparisons (`pairs`), the minimum significant
difference (`msd`, for equal replication) and the error used.

## Details

Methods: `"tukey"` (Tukey-Kramer HSD), `"lsd"` (Fisher), `"bonferroni"`,
`"sidak"`, `"holm"`, `"duncan"`, `"snk"` (Student-Newman-Keuls),
`"regwq"`, `"scheffe"`, `"dunnett"` (each level against `control`; the
groups are then `"*"` for a difference with the control and `"ns"`) and
`"games_howell"` (unequal variances; own SE and Welch df per pair).

With `by`, the comparison is made within each level of another factor
(simple effects). In split-plots the sub-plot factor within a main plot
uses Error b, and the main-plot factor within a sub-plot level uses the
pooled error `[(b - 1) Eb + Ea] / b` with Satterthwaite degrees of
freedom; strip-plots pool their errors the same way.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
compare_means(fit, "Nitrogen_kg_ha")
#> 
#> Tukey HSD (Tukey-Kramer), alpha = 0.05
#> Response: Yield_t_ha   compared: Nitrogen_kg_ha
#> Error: Residual error (MS = 0.01643, df = 9)
#> Minimum significant difference = 0.2829
#> 
#>  Nitrogen_kg_ha n  mean      se raw_mean     sd group
#>             180 4 7.183 0.06409    7.183 0.2943     a
#>             120 4 6.885 0.06409    6.885 0.3663     b
#>              60 4 5.845 0.06409    5.845 0.2978     c
#>               0 4 4.215 0.06409    4.215 0.3101     d
#> 
#> Means sharing a letter do not differ.

fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
compare_means(fit, "Variety", method = "dunnett", control = "Negro Jamapa")
#> 
#> Dunnett (vs control), alpha = 0.05
#> Response: Yield_kg_plot   compared: Variety
#> Error: Residual error (MS = 0.01054, df = 15)
#> Two-sided comparisons with the control (Negro Jamapa). 
#> Minimum significant difference = 0.198
#> 
#>         Variety n  mean      se raw_mean      sd   group
#>        Azufrado 4 2.812 0.05133    2.812 0.10626       *
#>    Flor de Mayo 4 2.597 0.05133    2.598 0.09708       *
#>  Pinto Saltillo 4 2.420 0.05133    2.420 0.09345       *
#>    Negro Jamapa 4 2.132 0.05133    2.132 0.10210 control
#>            Bayo 4 1.917 0.05133    1.917 0.11325       *

fit <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "split_rcbd")
compare_means(fit, "Variety", by = "Irrigation")
#> 
#> Tukey HSD (Tukey-Kramer), alpha = 0.05
#> Response: Yield_t_ha   compared: Variety within each level of Irrigation
#> Error: Error b (sub-plot) (MS = 0.0004931, df = 12)
#> Minimum significant difference = 0.05383
#> 
#>  Irrigation Variety n  mean      se raw_mean     sd group
#>        Drip      V2 3 9.123 0.01282    9.123 0.1716     a
#>        Drip      V4 3 8.913 0.01282    8.913 0.1818     b
#>        Drip      V1 3 8.387 0.01282    8.387 0.1823     c
#>        Drip      V3 3 7.810 0.01282    7.810 0.1706     d
#>      Furrow      V2 3 7.900 0.01282    7.900 0.1411     a
#>      Furrow      V4 3 7.617 0.01282    7.617 0.1861     b
#>      Furrow      V1 3 7.277 0.01282    7.277 0.1823     c
#>      Furrow      V3 3 6.853 0.01282    6.853 0.1320     d
#> 
#> Means sharing a letter do not differ (letters valid within each level).
```
