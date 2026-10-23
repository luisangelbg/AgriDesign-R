# Field book of a randomised layout

The field book of
[`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)
with empty columns for the responses to be measured, ready to be filled
in and analysed with
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

## Usage

``` r
field_book(layout, responses = "Response", file = NULL)
```

## Arguments

- layout:

  An `"ad_layout"` object.

- responses:

  Names of the response columns to add.

- file:

  Optional CSV file to write.

## Value

The field book as a data frame (invisibly when `file` is given).

## Examples

``` r
lay <- make_design("rcbd", trt = 4, r = 3, seed = 7)
head(field_book(lay, c("Yield", "Height")))
#>   Plot Block FieldRow FieldCol Treatment Yield Height
#> 1    1     1        1        1        T2    NA     NA
#> 2    2     1        1        2        T3    NA     NA
#> 3    3     1        1        3        T1    NA     NA
#> 4    4     1        1        4        T4    NA     NA
#> 5    5     2        2        4        T4    NA     NA
#> 6    6     2        2        3        T1    NA     NA
```
