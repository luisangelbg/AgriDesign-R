# Box-Cox profile likelihood of a design

Profile log-likelihood of the Box-Cox power `lambda` for the model of
the design, its maximum, the 95 % interval and the nearest conventional
power (-1, -0.5, 0, 0.5, 1 or 2).

## Usage

``` r
box_cox(fit, lambda = seq(-2, 2, by = 0.1))
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- lambda:

  Grid of powers.

## Value

A list with the profile (`profile`), the estimate (`lambda`), the 95 %
interval (`lower`, `upper`), the nearest conventional power (`rounded`)
and whether the interval contains 1 (`contains_one`); `NULL` when the
response has zeros or negative values.

## Examples

``` r
fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
                    block = "Block")
box_cox(fit)[c("lambda", "lower", "upper", "rounded")]
#> $lambda
#> [1] 0.2
#> 
#> $lower
#> [1] 0.1
#> 
#> $upper
#> [1] 0.4
#> 
#> $rounded
#> [1] 0
#> 
```
