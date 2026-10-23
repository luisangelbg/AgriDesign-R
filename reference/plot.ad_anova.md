# Residual plots of a design

Six diagnostic plots of the residuals: residuals against fitted values,
normal Q-Q plot of the studentized residuals with a 95 % pointwise band,
scale-location, residuals by treatment, residuals in data order and a
histogram with the normal curve. Observations with \|studentized
residual\| above 2.5 are labelled with their row number.

## Usage

``` r
# S3 method for class 'ad_anova'
plot(x, which = 1:6, ...)
```

## Arguments

- x:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- which:

  Plots to draw (1 to 6). Several plots share one page.

- ...:

  Not used.

## Value

Invisibly, a data frame with the fitted values, residuals and
studentized residuals.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
plot(fit)

plot(fit, which = 2)
```
