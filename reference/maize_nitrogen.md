# Maize yield under four nitrogen rates (RCBD)

Randomised complete block design with four nitrogen rates in four
blocks.

## Usage

``` r
maize_nitrogen
```

## Format

A data frame with 16 rows and 3 columns:

- Block:

  Block (1-4).

- Nitrogen_kg_ha:

  Nitrogen rate, kg/ha (0, 60, 120, 180).

- Yield_t_ha:

  Grain yield, t/ha.

## Examples

``` r
design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#> 
#> Randomised complete block design (RCBD)
#> Response: Yield_t_ha   (n = 16)
#> Type III sums of squares
#> 
#>                Df  Sum Sq Mean Sq F value    Pr(>F)    
#> Block           3  1.0691 0.35637   21.69 1.864e-04 ***
#> Nitrogen_kg_ha  3 21.5509 7.18362  437.27 4.591e-10 ***
#> Residuals       9  0.1479 0.01643                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Grand mean = 6.03188   CV = 2.12 %   R-squared = 0.9935
```
