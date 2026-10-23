# Least-squares means of a design

Means of the levels of one or more factors, adjusted for the other terms
of the model: the other factors get equal weights and covariates are
held at their mean.

## Usage

``` r
ls_means(fit, factors)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- factors:

  Name(s) of the factor(s) whose means are wanted.

## Value

A data frame with one row per level (or combination of levels): the
least-squares mean, its standard error, the number of observations and
the raw mean and standard deviation.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha",
                    block = "Block", design = "rcbd")
ls_means(fit, "Nitrogen_kg_ha")
#>   Nitrogen_kg_ha   mean        se n raw_mean        sd label
#> 1              0 4.2150 0.0640868 4   4.2150 0.3101075     0
#> 2             60 5.8450 0.0640868 4   5.8450 0.2978255    60
#> 3            120 6.8850 0.0640868 4   6.8850 0.3662877   120
#> 4            180 7.1825 0.0640868 4   7.1825 0.2943213   180
```
