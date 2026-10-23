# Standard methods for fitted designs

[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html) and
[`anova()`](https://rdrr.io/r/stats/anova.html) work on the result of
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
as on a linear model. [`anova()`](https://rdrr.io/r/stats/anova.html)
returns the analysis of variance as a data frame of class `"anova"`,
with the error stratum against which each term is tested.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
result of
[`compare_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/compare_means.md)
returns its table of means with their letters.

## Usage

``` r
# S3 method for class 'ad_anova'
residuals(object, type = c("response", "studentized"), ...)

# S3 method for class 'ad_anova'
fitted(object, ...)

# S3 method for class 'ad_anova'
nobs(object, ...)

# S3 method for class 'ad_anova'
anova(object, ...)

# S3 method for class 'ad_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- object, x:

  A fitted design from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
  (or, for
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), a
  comparison from
  [`compare_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/compare_means.md)).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (observed minus fitted) or `"studentized"` (internally
  studentized, as
  [`stats::rstandard()`](https://rdrr.io/r/stats/influence.measures.html)).

- row.names, optional, ...:

  Passed on or not used.

## Value

A numeric vector (residuals, fitted values), a count (`nobs`), or a data
frame.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
anova(fit)
#> Analysis of variance: Randomised complete block design (RCBD)
#> Response: Yield_t_ha
#> 
#>                Df  Sum Sq Mean Sq F value     Pr(>F) Tested against
#> Block           3  1.0691  0.3564  21.692 0.00018644              2
#> Nitrogen_kg_ha  3 21.5509  7.1836 437.267 0.00000000              2
#> Residuals       9  0.1479  0.0164                                 1
head(residuals(fit, type = "studentized"))
#>           1           2           3           4           5           6 
#>  0.25356235 -0.16253997  0.77369025 -0.86471264  0.14953677 -0.05851439 
nobs(fit)
#> [1] 16
as.data.frame(compare_means(fit, "Nitrogen_kg_ha"))
#>   Nitrogen_kg_ha   mean        se n raw_mean        sd label group
#> 1              0 4.2150 0.0640868 4   4.2150 0.3101075     0     d
#> 2             60 5.8450 0.0640868 4   5.8450 0.2978255    60     c
#> 3            120 6.8850 0.0640868 4   6.8850 0.3662877   120     b
#> 4            180 7.1825 0.0640868 4   7.1825 0.2943213   180     a
```
