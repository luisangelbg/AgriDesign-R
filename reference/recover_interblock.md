# Recover inter-block information (REML)

In incomplete block designs (alpha, lattice, BIBD, cyclic, augmented),
the intra-block analysis of
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
treats blocks as fixed and uses only comparisons within blocks. When
blocks are random, block totals also carry information on the
treatments. This function fits blocks as a random effect, estimates the
block and error variances by restricted maximum likelihood (REML) and
returns the treatment means as generalised least-squares estimates
(BLUEs), with their standard errors and the efficiency gained over the
intra-block analysis and over a complete block design.

## Usage

``` r
recover_interblock(fit)
```

## Arguments

- fit:

  An `"ad_anova"` object of an incomplete block design (designs `"ibd"`,
  `"bib"`, `"augmented"`, `"rcbd"` with incomplete blocks).

## Value

An object of class `"ad_recovery"` with the variance components
(`sigma2_block`, `sigma2`), the adjusted treatment means and SEs
(`means`), the average standard error of a difference for the combined,
intra-block and complete-block analyses, and the REML log-likelihood.

## Examples

``` r
lay <- make_design("alpha", trt = 12, r = 3, k = 3, seed = 4)
b <- field_book(lay, "Yield")
b$Yield <- 10 + as.integer(factor(b$Treatment)) / 5 + rep(c(-1, 0, 1), length.out = nrow(b))
fit <- design_anova(b, "Yield", "Treatment", block = c("Block", "IncBlock"))
recover_interblock(fit)
#> 
#> Recovery of inter-block information (REML)
#> Incomplete block design (lattice, alpha, BIBD), intra-block analysis
#> Variance components: blocks = 0, error = 0.7576 (block variance on the boundary: no inter-block information)
#> Average SE of a difference: combined 0.7107, intra-block 1.063, replicates as complete blocks 0.7107
#> Relative efficiency of the combined analysis: 223.9 % of the intra-block analysis, 100.0 % of complete blocks
#> 
#>  level n  mean     se intra_mean
#>    T11 3 10.27 0.5025     10.287
#>     T4 3 11.40 0.5025     11.420
#>     T6 3 12.13 0.5025     12.211
#>     T5 3 10.93 0.5025     10.763
#>    T10 3 10.73 0.5025     10.730
#>     T1 3 10.87 0.5025     11.039
#>     T3 3 11.53 0.5025     11.573
#>     T9 3 12.07 0.5025     12.249
#>    T12 3 11.47 0.5025     11.639
#>     T7 3 11.67 0.5025     11.635
#>     T2 3 10.33 0.5025      9.911
#>     T8 3 12.20 0.5025     12.144
```
