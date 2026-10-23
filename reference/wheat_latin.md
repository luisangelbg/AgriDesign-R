# Wheat yield in a 4 x 4 Latin square

Wheat yield in a 4 x 4 Latin square

## Usage

``` r
wheat_latin
```

## Format

A data frame with 16 rows and 4 columns:

- Row:

  Row of the field (1-4).

- Column:

  Column of the field (1-4).

- Treatment:

  Treatment (A-D).

- Yield_t_ha:

  Grain yield, t/ha.

## Examples

``` r
design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column")
#> 
#> Latin square
#> Response: Yield_t_ha   (n = 16)
#> Type III sums of squares
#> 
#>           Df  Sum Sq  Mean Sq F value    Pr(>F)    
#> Row        3 0.05885 0.019617   2.812 1.301e-01    
#> Column     3 0.02165 0.007217   1.035 4.420e-01    
#> Treatment  3 3.43375 1.144583 164.098 3.802e-06 ***
#> Residuals  6 0.04185 0.006975                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Grand mean = 3.9075   CV = 2.14 %   R-squared = 0.9882
```
