# Yield of five bean varieties (CRD)

Completely randomised design with five common bean varieties.

## Usage

``` r
bean_varieties
```

## Format

A data frame with 20 rows and 3 columns:

- Plot:

  Plot number.

- Variety:

  Bean variety.

- Yield_kg_plot:

  Yield, kg per plot.

## Examples

``` r
design_anova(bean_varieties, "Yield_kg_plot", "Variety")
#> 
#> Completely randomised design (CRD)
#> Response: Yield_kg_plot   (n = 20)
#> Type III sums of squares
#> 
#>           Df Sum Sq Mean Sq F value    Pr(>F)    
#> Variety    4 2.0442 0.51105   48.49 2.096e-08 ***
#> Residuals 15 0.1581 0.01054                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Grand mean = 2.376   CV = 4.32 %   R-squared = 0.9282
```
