# Sorghum hybrids in a balanced incomplete block design (simulated)

Simulated grain yields of seven sorghum hybrids in seven blocks of three
plots: every hybrid appears three times and every pair of hybrids meets
once in a block (t = 7, k = 3, r = 3, lambda = 1).

## Usage

``` r
bib_sorghum
```

## Format

A data frame with 21 rows and 3 columns:

- Block:

  Block (1-7).

- Hybrid:

  Hybrid (H1-H7).

- Yield_t_ha:

  Grain yield, t/ha.

## Examples

``` r
durbin_test(bib_sorghum, "Yield_t_ha", "Hybrid", "Block")
#> 
#> Durbin rank test for balanced incomplete blocks
#> Response: Yield_t_ha   factor: Hybrid
#> chisq = 10.44, df = 6, p = 0.1071
#> 
#> Pairwise comparisons: Conover t on rank sums (p adjusted by holm)
#>  level_i level_j statistic        p  p_adj different
#>       H1      H2   -0.3118 0.763160 1.0000     FALSE
#>       H1      H4   -1.8708 0.098282 1.0000     FALSE
#>       H1      H3    1.2472 0.247591 1.0000     FALSE
#>       H1      H5   -0.3118 0.763160 1.0000     FALSE
#>       H1      H6    1.8708 0.098282 1.0000     FALSE
#>       H1      H7   -0.6236 0.550247 1.0000     FALSE
#>       H2      H4   -1.5590 0.157609 1.0000     FALSE
#>       H2      H3    1.5590 0.157609 1.0000     FALSE
#>       H2      H5    0.0000 1.000000 1.0000     FALSE
#>       H2      H6    2.1826 0.060613 1.0000     FALSE
#>       H2      H7   -0.3118 0.763160 1.0000     FALSE
#>       H4      H3    3.1180 0.014271 0.2854     FALSE
#>       H4      H5    1.5590 0.157609 1.0000     FALSE
#>       H4      H6    3.7417 0.005692 0.1195     FALSE
#>       H4      H7    1.2472 0.247591 1.0000     FALSE
#>       H3      H5   -1.5590 0.157609 1.0000     FALSE
#>       H3      H6    0.6236 0.550247 1.0000     FALSE
#>       H3      H7   -1.8708 0.098282 1.0000     FALSE
#>       H5      H6    2.1826 0.060613 1.0000     FALSE
#>       H5      H7   -0.3118 0.763160 1.0000     FALSE
#>       H6      H7   -2.4944 0.037264 0.7080     FALSE
#> 
#> Groups (levels sharing a letter do not differ at alpha = 0.05):
#>  level n median  mean mean_rank rank_sum group
#>     H4 3   7.97 7.763     3.000      9.0     a
#>     H7 3   7.41 7.527     2.333      7.0     a
#>     H2 3   6.45 6.507     2.167      6.5     a
#>     H5 3   5.68 6.417     2.167      6.5     a
#>     H1 3   6.78 6.497     2.000      6.0     a
#>     H3 3   6.43 6.127     1.333      4.0     a
#>     H6 3   5.98 5.707     1.000      3.0     a
#> 
#> Design: 7 treatments in 7 blocks of 3 plots, r = 3, lambda = 1
#> F form: F = 3.917 on 6 and 8 df, p = 0.03973
```
