# Analysis of variance of an experimental design

Fits the linear model of a classical agricultural design and tests every
term against the right error stratum (for instance, main plots of a
split-plot against Error a).

## Usage

``` r
design_anova(
  data,
  response,
  treatment,
  block = NULL,
  row = NULL,
  col = NULL,
  covariate = NULL,
  design = NULL,
  ss_type = 3,
  random_blocks = TRUE,
  random_nested = TRUE,
  transform = NULL
)
```

## Arguments

- data:

  A data frame in long format: one row per plot.

- response:

  Name of the numeric response column.

- treatment:

  Name(s) of the treatment factor(s). In split-plots, strip-plots and
  nested designs the order matters: the first factor goes on main plots
  (or is the outer factor).

- block:

  Name of the block column; two names (replicate, block within
  replicate) for incomplete block designs.

- row, col:

  Names of the row and column factors of a Latin square.

- covariate:

  Name(s) of numeric covariates (analysis of covariance).

- design:

  Identifier of the design (see
  [`designs()`](https://luisangelbg.github.io/AgriDesign-R/reference/designs.md)).
  When `NULL` the design is suggested from the roles given.

- ss_type:

  Type of sums of squares: `3` (default, adjusted for all other terms)
  or `1` (sequential).

- random_blocks:

  In a generalised RCBD, treat blocks as random and test treatments
  against the block x treatment interaction.

- random_nested:

  In a nested design, treat the inner factor as random and test the
  outer factor against it.

- transform:

  Optional transformation of the response before the analysis, by its
  identifier in
  [`transformations()`](https://luisangelbg.github.io/AgriDesign-R/reference/transformations.md)
  (`"log"`, `"sqrt"`, `"asin"`...). Means are then also reported
  back-transformed.

## Value

An object of class `"ad_anova"`: a list with the ANOVA table (`table`),
the error strata (`strata`), the grand mean, the coefficient of
variation, R squared, the fitted model and the data used.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha",
                    block = "Block", design = "rcbd")
fit
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
