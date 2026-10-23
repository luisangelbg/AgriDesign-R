# Compare transformations of the response

Refits the design after each transformation that suits the range of the
response and repeats the assumption checks. A transformation is
suggested only when the original scale fails some check and the
transformed scale passes all of them: the conventional power nearest to
the Box-Cox estimate when it passes, otherwise the candidate with the
least evidence against the assumptions.

## Usage

``` r
transformations(fit, alpha = 0.05)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- alpha:

  Significance level of the checks.

## Value

An object of class `"ad_transformations"`: the table of candidates
(`table`), the identifier of the suggested transformation (`suggested`,
`NULL` when none is needed or none works), a short `verdict` and the
Box-Cox profile (`box_cox`).

## Examples

``` r
fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
                    block = "Block")
transformations(fit)
#> 
#> Transformations of Yield_kg_plant (same model, alpha = 0.05)
#> 
#>     Transformation          W     p norm. p Levene p addit. Fmax Skew  CV %
#>     None (original scale)   0.959 0.36    0.85     <1e-04   4.00  0.50 0.8 
#>  *  Natural logarithm ln(y) 0.958 0.34    0.94     0.220    1.85 -0.47 0.5 
#>     ln(y + 1)               0.968 0.56    0.52     0.535    1.78  0.16 0.3 
#>     log10(y)                0.958 0.34    0.94     0.220    1.85 -0.47 0.5 
#>     Square root sqrt(y)     0.946 0.17    0.77     0.018    2.02  0.73 0.3 
#>     sqrt(y + 0.5)           0.938 0.11    0.84     0.005    2.16  0.77 0.2 
#>     Arcsine sqrt(y/100)     0.945 0.16    0.78     0.013    2.05  0.74 0.3 
#>     Reciprocal 1/y          0.989 0.99    0.83     <1e-04   6.06 -0.01 1.3 
#> (p norm.: Shapiro-Wilk; p Levene: Levene with median; p addit.: Tukey's additivity test)
#> 
#> Box-Cox: lambda = 0.20, 95 % interval 0.1 to 0.4; nearest conventional power 0.
#> 
#> Suggested: Natural logarithm ln(y). It is the power closest to the Box-Cox estimate and it satisfies the checks (smallest p = 0.2204 vs 2.346e-06 on the original scale). Means are compared on the transformed scale; report them back-transformed.
```
