# Check the assumptions of the analysis of variance

Runs the usual checks on the residuals of a fitted design: normality
(Shapiro-Wilk, Anderson-Darling, Jarque-Bera), homogeneity of variances
among treatments (Levene with median and with mean, Bartlett,
Fligner-Killeen, Hartley's Fmax), additivity of blocks and treatments
(Tukey's one-degree-of-freedom test, only for blocked designs),
independence in data order (Durbin-Watson) and extreme studentized
residuals. The checks are summarised in a grade from A (assumptions
satisfied) to D (several assumptions fail).

## Usage

``` r
check_assumptions(fit, alpha = 0.05)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- alpha:

  Significance level of the checks.

## Value

An object of class `"ad_assumptions"` with the table of tests (`tests`),
the residuals and studentized residuals (`residuals`), the observations
with \|studentized residual\| \> 2.5 (`suspicious`) and the grade.

## Details

Each check is rated "ok" when its p-value is at least `alpha`, "warn"
when it lies between `alpha / 5` and `alpha` and "bad" below
`alpha / 5`. The grade is D with two or more "bad" checks, C with one
"bad" or two "warn", B with a single "warn" and A otherwise.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
check_assumptions(fit)
#> 
#> Assumptions of the ANOVA for Yield_t_ha (alpha = 0.05)
#> 
#>  Check        Test                        Statistic    df    p      Rating
#>  normality    Shapiro-Wilk                W = 0.9304         0.2477 ok    
#>  normality    Anderson-Darling            A2 = 0.4947        0.1840       
#>  normality    Jarque-Bera                 JB = 2.303   2     0.3162       
#>  variances    Levene (median)             F = 0.5348   3, 12 0.6672 ok    
#>  variances    Levene (mean)               F = 4.27     3, 12 0.0287       
#>  variances    Bartlett                    K2 = 0.1708  3     0.9822       
#>  variances    Fligner-Killeen             X2 = 1.314   3     0.7259       
#>  variances    Hartley Fmax                Fmax = 1.549                    
#>  additivity   Tukey one-degree-of-freedom F = 0.0548   1, 8  0.8208 ok    
#>  extremes     |studentized residual| > 3  count = 0                 ok    
#>  independence Durbin-Watson (data order)  d = 2.451                 ok    
#> 
#> Observations with |studentized residual| > 2.5:
#>  row group observed fitted studentized
#>   16   180     7.31  7.056       2.646
#> 
#> Grade A. Assumptions satisfied. Proceed with the analysis of variance on this scale.
```
