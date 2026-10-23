# Partition of the total sum of squares

A bar split into the share of the total sum of squares taken by every
source of the analysis of variance, including the error strata and the
residual.

## Usage

``` r
plot_ss_partition(fit, col = NULL, main = NULL)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- col:

  Colours of the sources.

- main:

  Title.

## Value

Invisibly, a data frame with the sum of squares and percentage of every
source.

## Examples

``` r
fit <- design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
                    block = "Block", design = "split_rcbd")
plot_ss_partition(fit)
```
