# Line x tester trial in maize

Simulated yields of six lines crossed with three testers, evaluated with
the nine parents in three randomised complete blocks. Parents have an
empty tester (lines) or an empty line (testers).

## Usage

``` r
maize_lxt
```

## Format

A data frame with 81 rows and 4 columns:

- Block:

  Block (1-3).

- Line:

  Line (L1-L6), `NA` for the tester parents.

- Tester:

  Tester (T1-T3), `NA` for the line parents.

- Yield_t_ha:

  Grain yield, t/ha.

## Examples

``` r
line_tester(maize_lxt, "Yield_t_ha", "Line", "Tester", block = "Block")
#> 
#> Line x tester analysis (Kempthorne 1957)
#> Response: Yield_t_ha   6 lines x 3 testers, 9 parents, r = 3
#> 
#>                       Df   Sum Sq   Mean Sq  F value    Pr(>F)    
#> Blocks                 2   0.6015   0.30076    5.880 4.988e-03  **
#> Entries               26 128.2189   4.93150   96.409 2.175e-35 ***
#>   Parents              8   3.8639   0.48299    9.442 6.072e-08 ***
#>     Lines (parents)    5   3.0794   0.61588   12.040 9.132e-08 ***
#>     Testers (parents)  2   0.5283   0.26413    5.164 9.004e-03  **
#>     Lines vs testers   1   0.2563   0.25627    5.010 2.951e-02   *
#>   Parents vs crosses   1 107.0510 107.05096 2092.815 1.114e-43 ***
#>   Crosses             17  17.3040   1.01788   19.899 5.536e-17 ***
#>     Lines              5  10.0670   2.01340   39.361 1.609e-16 ***
#>     Testers            2   2.8244   1.41219   27.608 6.751e-09 ***
#>     Lines x Testers   10   4.4126   0.44126    8.626 4.125e-08 ***
#> Error                 52   2.6599   0.05115                       
#> Total                 80 131.4803   1.64350                       
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Line x tester analysis (random model: lines and testers against L x T):
#>                 Df Sum Sq Mean Sq F value    Pr(>F)    
#> Lines            5 10.067 2.01340   4.563 1.989e-02   *
#> Testers          2  2.824 1.41219   3.200 8.427e-02   .
#> Lines x Testers 10  4.413 0.44126   8.626 4.125e-08 ***
#> Error           52  2.660 0.05115                      
#> 
#> Tested against another mean square:
#>   Lines vs Lines x Testers
#>   Testers vs Lines x Testers
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> GCA effects of lines:
#>  line  mean     gca      se      t         p
#>    L1 7.512  0.4913 0.07539  6.517 2.915e-08
#>    L2 6.691 -0.3298 0.07539 -4.375 5.868e-05
#>    L3 7.207  0.1857 0.07539  2.464 1.709e-02
#>    L4 6.321 -0.6998 0.07539 -9.283 1.274e-12
#>    L5 7.499  0.4780 0.07539  6.340 5.571e-08
#>    L6 6.896 -0.1254 0.07539 -1.663 1.023e-01
#> 
#> GCA effects of testers:
#>  tester  mean     gca      se     t         p
#>      T1 7.333  0.3124 0.05331  5.86 3.196e-07
#>      T2 6.792 -0.2287 0.05331 -4.29 7.781e-05
#>      T3 6.937 -0.0837 0.05331 -1.57 1.224e-01
#> 
#> SCA effects:
#>       T1    T2    T3
#> L1  0.19 -0.27  0.08
#> L2  0.03 -0.24  0.21
#> L3 -0.04  0.44 -0.40
#> L4 -0.42  0.54 -0.11
#> L5 -0.14 -0.21  0.36
#> L6  0.39 -0.25 -0.14
#> 
#> Standard errors:
#>              GCA line            GCA tester                   SCA 
#>               0.07539               0.05331               0.13060 
#>   GCA line difference GCA tester difference        SCA difference 
#>               0.10660               0.07539               0.18470 
#> 
#> Genetic components:
#>              component estimate
#>         Cov HS (lines)  0.17468
#>       Cov HS (testers)  0.05394
#>       Cov HS (average)  0.12102
#>                 Cov FS  0.29158
#>  sigma2 GCA (= Cov HS)  0.12102
#>             sigma2 SCA  0.13004
#> 
#>  F sigma2_A sigma2_D dominance
#>  0   0.4841   0.5201     1.466
#>  1   0.2420   0.1300     1.037
#> 
#> Contribution to the variation among crosses (%):
#>   lines testers     lxt 
#>   58.18   16.32   25.50 
```
