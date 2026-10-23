# Tomato varieties by fertiliser (factorial in RCBD)

Three tomato varieties by three fertiliser levels in three blocks.

## Usage

``` r
tomato_factorial
```

## Format

A data frame with 27 rows and 5 columns:

- Block:

  Block (1-3).

- Variety:

  Tomato variety.

- Fertilizer:

  Fertiliser level (F0, F1, F2).

- Yield_kg_plant:

  Fruit yield, kg per plant.

- Height_cm:

  Plant height, cm.

## Examples

``` r
design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
             block = "Block")
#> 
#> Factorial in RCBD
#> Response: Yield_kg_plant   (n = 27)
#> Type III sums of squares
#> 
#>                    Df   Sum Sq   Mean Sq F value    Pr(>F)    
#> Block               2  0.21509 0.1075444   152.1 3.882e-11 ***
#> Variety             2 12.29736 6.1486778  8697.5 5.086e-25 ***
#> Fertilizer          2  4.10702 2.0535111  2904.8 3.238e-21 ***
#> Variety:Fertilizer  4  0.49009 0.1225222   173.3 5.916e-13 ***
#> Residuals          16  0.01131 0.0007069                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Grand mean = 3.15556   CV = 0.84 %   R-squared = 0.9993
```
