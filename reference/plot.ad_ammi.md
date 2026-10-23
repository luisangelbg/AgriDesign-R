# Biplots of an AMMI analysis

The AMMI1 biplot puts the means of genotypes and environments on the
horizontal axis and their first interaction score on the vertical axis;
the AMMI2 biplot shows the first two interaction axes, with the
environments as vectors from the origin.

## Usage

``` r
# S3 method for class 'ad_ammi'
plot(
  x,
  type = c("ammi1", "ammi2"),
  col = NULL,
  labels = TRUE,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  cex = 0.8,
  ...
)
```

## Arguments

- x:

  An `"ad_ammi"` object from
  [`ammi()`](https://luisangelbg.github.io/AgriDesign-R/reference/ammi.md).

- type:

  `"ammi1"` (mean vs IPCA1) or `"ammi2"` (IPCA1 vs IPCA2).

- col:

  Two colours: genotypes and environments.

- labels:

  Write the names of genotypes and environments.

- main, xlab, ylab:

  Titles; `NULL` builds them from the analysis.

- cex:

  Size of the labels.

- ...:

  Further graphical parameters passed to
  [`graphics::plot.window()`](https://rdrr.io/r/graphics/plot.window.html).

## Value

Invisibly, a data frame with the coordinates of every point.

## Examples

``` r
fit <- ammi(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
plot(fit, type = "ammi2")
```
