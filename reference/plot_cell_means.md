# Heat map of cell means

Least-squares means of every combination of two treatment factors as a
coloured table with the values written in the cells.

## Usage

``` r
plot_cell_means(fit, rows, cols, digits = 2, main = NULL)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- rows, cols:

  Factors in the rows and columns of the map.

- digits:

  Decimals of the values.

- main:

  Title; `NULL` builds it from the analysis.

## Value

Invisibly, the cell means plotted.

## Examples

``` r
fit <- design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
                    block = "Block")
plot_cell_means(fit, "Variety", "Fertilizer")
```
