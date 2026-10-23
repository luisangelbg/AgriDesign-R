# Trend analysis of a quantitative factor

Splits the sum of squares of a factor with numeric levels (doses, rates,
dates) into linear, quadratic, cubic... components with orthogonal
polynomial contrasts (also with unequal spacing) and fits the polynomial
of the highest significant degree to the observations. A fitted
quadratic reports where the maximum (or minimum) response lies.

## Usage

``` r
trend(fit, factor, degree = NULL, alpha = 0.05)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- factor:

  Treatment factor whose means are compared; two names compare all their
  combinations together (not available for split- or strip-plots, where
  `by` should be used).

- degree:

  Highest degree to test (default: up to 4, at most the number of levels
  minus one).

- alpha:

  Significance level.

## Value

An object of class `"ad_trend"` with the table of components
(`components`), the fitted coefficients (`coef`, intercept first), its R
squared and the optimum of a quadratic (`optimum`).

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
trend(fit, "Nitrogen_kg_ha")
#> 
#> Trend analysis of Yield_t_ha on Nitrogen_kg_ha (error: Residual error)
#> 
#>  component coefficients        ss         F df_error         p
#>     Linear -3, -1, 1, 3 19.770661 1203.4388        9 6.788e-11
#>  Quadratic 1, -1, -1, 1  1.775556  108.0780        9 2.587e-06
#>      Cubic -1, 3, -3, 1  0.004651    0.2831        9 6.075e-01
#> 
#> Highest significant component: quadratic.
#> Fitted response: y = 4.207 + 0.03323 x - 9.253e-05 x^2   R-squared = 0.946
#> Maximum response at x = 179.5 (y = 7.19); at the edge of the levels tested
```
