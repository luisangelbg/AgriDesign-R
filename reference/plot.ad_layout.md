# Plot a randomised field layout

Draws the plots on their field rows and columns coloured by treatment
(or by main-plot factor), with the plot numbers and the outlines of
blocks and main plots.

## Usage

``` r
# S3 method for class 'ad_layout'
plot(x, color_by = NULL, labels = TRUE, numbers = TRUE, main = NULL, ...)
```

## Arguments

- x:

  An `"ad_layout"` object from
  [`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md).

- color_by:

  Column of the field book that colours the plots (default: `A` in split
  and factorial designs, otherwise `Treatment`).

- labels:

  Write the treatment in each plot.

- numbers:

  Write the plot number in each plot.

- main:

  Title.

- ...:

  Not used.

## Value

Invisibly, the field book.

## Examples

``` r
plot(make_design("strip", A = c("Drip", "Furrow", "Sprinkler"), B = c("V1", "V2"), r = 2, seed = 3))
```
