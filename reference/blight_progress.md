# Late blight progress under four fungicide programmes (simulated)

Simulated severity of late blight (percentage of leaf area affected) in
a potato trial with four fungicide programmes in four blocks, assessed
weekly from day 21 to day 63 after planting. The values were generated
with logistic epidemics (data-raw/make-data-pathology.R).

## Usage

``` r
blight_progress
```

## Format

A data frame with 112 rows and 4 columns:

- Block:

  Block (1-4).

- Treatment:

  Fungicide programme: Untreated, Contact, Systemic, Alternated.

- Day:

  Days after planting.

- Severity_pct:

  Disease severity, percent.

## Examples

``` r
a <- audpc_by_unit(blight_progress, "Severity_pct", "Day", c("Block", "Treatment"))
design_anova(a, "audpc", "Treatment", block = "Block")
#> 
#> Randomised complete block design (RCBD)
#> Response: audpc   (n = 16)
#> Type III sums of squares
#> 
#>           Df  Sum Sq Mean Sq F value    Pr(>F)    
#> Block      3  210707   70236    4.16 4.179e-02   *
#> Treatment  3 7057193 2352398  139.33 7.348e-08 ***
#> Residuals  9  151948   16883                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Grand mean = 2247.53   CV = 5.78 %   R-squared = 0.9795
```
