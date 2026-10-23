# Irrigation systems by varieties (split-plot in RCBD)

Two irrigation systems on main plots and four varieties on sub-plots, in
three blocks.

## Usage

``` r
split_irrigation
```

## Format

A data frame with 24 rows and 4 columns:

- Block:

  Block (1-3).

- Irrigation:

  Irrigation system (main-plot factor).

- Variety:

  Variety (sub-plot factor).

- Yield_t_ha:

  Yield, t/ha.

## Examples

``` r
design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
             block = "Block", design = "split_rcbd")
#> 
#> Split-plot in RCBD (main plots in blocks)
#> Response: Yield_t_ha   (n = 24)
#> Type III sums of squares
#> 
#>                     Df   Sum Sq   Mean Sq  F value    Pr(>F)    
#> Block                2 0.452575 0.2262875   300.05 3.322e-03  **
#> Irrigation           1 7.889067 7.8890667 10460.64 9.558e-05 ***
#> Error a (main plot)  2 0.001508 0.0007542                       
#> Variety              3 4.836800 1.6122667  3269.95 9.749e-18 ***
#> Irrigation:Variety   3 0.098733 0.0329111    66.75 9.343e-08 ***
#> Error (residual)    12 0.005917 0.0004931                       
#> 
#> Tested against a higher error stratum:
#>   Block vs Error a (main plot)
#>   Irrigation vs Error a (main plot)
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Grand mean = 7.985   CV = 0.28 %   R-squared = 0.9996
```
