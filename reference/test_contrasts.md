# Test contrasts among the levels of a factor

A contrast compares groups of treatments with coefficients that sum to
zero, one per level in the order of `levels(fit$data[[factor]])`; for
instance a control against the other three treatments is
`c(3, -1, -1, -1)`.

## Usage

``` r
test_contrasts(fit, factor, coef)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- factor:

  Treatment factor whose means are compared; two names compare all their
  combinations together (not available for split- or strip-plots, where
  `by` should be used).

- coef:

  A vector of coefficients, or a matrix with one contrast per row (row
  names are used as labels).

## Value

A data frame with the estimate, its SE, t, p, the sum of squares of the
contrast and F (1 and `df_error` degrees of freedom).

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
test_contrasts(fit, "Nitrogen_kg_ha", rbind("0 vs fertilised" = c(3, -1, -1, -1),
                                      "120 vs 180" = c(0, 0, 1, -1)))
#>          contrast estimate         se         t            p         ss
#> 1 0 vs fertilised  -7.2675 0.22200319 -32.73602 1.139670e-10 17.6055187
#> 2      120 vs 180  -0.2975 0.09063242  -3.28249 9.491459e-03  0.1770125
#>            F df_error
#> 1 1071.64674        9
#> 2   10.77474        9
```
