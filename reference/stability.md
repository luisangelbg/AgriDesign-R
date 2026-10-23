# Stability statistics of a multi-environment trial

Parametric and rank-based stability statistics computed from the table
of genotype means in every environment and the pooled error of the
trial:

- `"eberhart"`:

  Eberhart and Russell (1966): regression coefficient `b` on the
  environmental index (with its standard error from the genotype's own
  deviations, and a t test of `b = 1`), deviation mean square `s2d`
  (mean scale, minus the error of a mean) with its F test against the
  pooled error, and their analysis of variance.

- `"finlay"`:

  Finlay and Wilkinson (1963): regression on the environment mean,
  optionally on the log scale.

- `"shukla"`:

  Shukla's (1972) stability variance `sigma2` and the same after
  removing the linear regression on the environmental index (`s2`), with
  F tests against the pooled error (plot scale).

- `"wricke"`:

  Wricke's (1962) ecovalence (plot scale) and its share of the
  interaction.

- `"kang"`:

  Kang's (1993) yield-stability statistic YS: the rank of the mean
  (highest = g), adjusted by the LSD, plus the stability rating from
  Shukla's test (-8, -4, -2 at 0.01, 0.05, 0.10).

- `"huehn"`:

  Huehn's (1979) rank statistics S1 and S2 on the ranks corrected for
  the genotype effects, with the Z tests of Nassar and Huehn (1987).

The common analysis of variance follows Shukla (1972): the interaction
is split into the heterogeneity of the regressions and the residual,
with (g - 1)(e - 2) degrees of freedom (Eberhart and Russell printed
g(e - 2), but e - 2 of them belong to environments and are identically
zero).

## Usage

``` r
stability(
  data,
  response,
  genotype,
  environment,
  block = NULL,
  mse = NULL,
  df_error = NULL,
  reps = 1,
  method = c("eberhart", "finlay", "shukla", "wricke", "kang", "huehn"),
  alpha = 0.05,
  log_finlay = FALSE
)
```

## Arguments

- data:

  A data frame in long format.

- response:

  Name of the numeric response column.

- genotype, environment:

  Names of the genotype and environment columns.

- block:

  Optional name of the block (replicate) column; blocks are taken as
  nested within environments.

- mse, df_error:

  Error mean square (plot scale) and its degrees of freedom, only when
  the data hold one mean per cell.

- reps:

  Number of replicates behind each mean, when the data hold one mean per
  cell.

- method:

  One or more of `"eberhart"`, `"finlay"`, `"shukla"`, `"wricke"`,
  `"kang"` and `"huehn"` (default: all).

- alpha:

  Significance level of the LSD in Kang's statistic and of the critical
  values of Huehn's tests.

- log_finlay:

  Fit the regression of Finlay and Wilkinson on the natural logarithm of
  the means, as in their paper.

## Value

An object of class `"ad_stability"`: the analysis of variance (`anova`)
and one element per method.

## References

Eberhart, S. A. and Russell, W. A. (1966) Crop Science 6, 36-40. Finlay,
K. W. and Wilkinson, G. N. (1963) Australian Journal of Agricultural
Research 14, 742-754. Shukla, G. K. (1972) Heredity 29, 237-245. Wricke,
G. (1962) Zeitschrift fur Pflanzenzuchtung 47, 92-96. Kang, M. S. (1993)
Agronomy Journal 85, 754-757. Huehn, M. (1979) EDV in Medizin und
Biologie 10, 112-117. Nassar, R. and Huehn, M. (1987) Biometrics 43,
45-53.

## Examples

``` r
st <- stability(maize_met, "Yield_t_ha", "Hybrid", "Environment", block = "Block")
st
#> 
#> Stability analysis of a multi-environment trial
#> Response: Yield_t_ha   12 genotypes x 6 environments, 3 replicate(s) per cell
#> 
#> Analysis of variance (plot scale; Shukla 1972):
#>                  Df  Sum Sq Mean Sq  F value    Pr(>F)    
#> Genotypes        11  55.358  5.0326   5.9765 2.647e-06 ***
#> Environments      5 223.530 44.7060 496.7562 9.569e-84 ***
#> GxE              55  46.314  0.8421   9.3567 7.883e-26 ***
#>   Heterogeneity  11   1.786  0.1624   0.1604 9.988e-01    
#>   Residual       44  44.528  1.0120  11.2449 1.286e-27 ***
#> Pooled error    132  11.879  0.0900                       
#> Total            71 325.202  4.5803                       
#> 
#> Tested against another mean square:
#>   Genotypes vs GxE
#>   Heterogeneity vs Residual
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Eberhart and Russell (1966): regression on the environmental index
#>                          Df  Sum Sq  Mean Sq F value    Pr(>F)    
#> Genotypes                11  55.358   5.0326  4.9729 5.596e-05 ***
#> Env + (G x E)            60 269.844   4.4974                      
#>   Environments (linear)   1 223.530 223.5299                      
#>   G x E (linear)         11   1.786   0.1624  0.1604 9.988e-01    
#>   Pooled deviation       44  44.528   1.0120 11.2449 1.286e-27 ***
#> Pooled error            132  11.879   0.0900                      
#> 
#> Tested against another mean square:
#>   Genotypes vs Pooled deviation
#>   G x E (linear) vs Pooled deviation
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  genotype  mean      b    se_b     t_b1    p_b1       s2d   F_dev     p_dev
#>       H01 6.744 1.0761 0.47231  0.16108 0.87984  1.355129 46.1731 3.253e-24
#>       H02 7.810 1.1125 0.17013  0.66155 0.54442  0.149721  5.9909 1.842e-04
#>       H03 7.151 1.0137 0.19387  0.07069 0.94704  0.203370  7.7793 1.164e-05
#>       H04 6.728 0.9249 0.08022 -0.93652 0.40203  0.009962  1.3321 2.614e-01
#>       H05 6.697 1.0482 0.25099  0.19199 0.85711  0.361168 13.0395 5.618e-09
#>       H06 7.547 1.0506 0.19287  0.26257 0.80585  0.200983  7.6997 1.314e-05
#>       H07 8.273 1.0509 0.17414  0.29224 0.78463  0.158303  6.2770 1.178e-04
#>       H08 6.832 0.9708 0.19796 -0.14760 0.88980  0.213339  8.1116 7.032e-06
#>       H09 7.846 1.1084 0.17337  0.62540 0.56562  0.156641  6.2216 1.285e-04
#>       H10 7.622 0.8438 0.07592 -2.05768 0.10874  0.005788  1.1930 3.169e-01
#>       H11 7.208 0.9597 0.30665 -0.13145 0.90177  0.553866 19.4630 1.312e-12
#>       H12 7.731 0.8404 0.04465 -3.57475 0.02328 -0.017620  0.4126 7.993e-01
#>      r2
#>  0.5648
#>  0.9145
#>  0.8724
#>  0.9708
#>  0.8134
#>  0.8812
#>  0.9010
#>  0.8574
#>  0.9109
#>  0.9686
#>  0.7100
#>  0.9888
#> b: SE from the genotype's own deviations (e - 2 df), t tests b = 1; s2d = MS deviations - error of a mean.
#> 
#> Finlay and Wilkinson (1963): regression on the environment mean  
#>  genotype  mean intercept      b    se_b     r2
#>       H01 6.744   -1.1638 1.0761 0.47231 0.5648
#>       H02 7.810   -0.3663 1.1125 0.17013 0.9145
#>       H03 7.151   -0.2987 1.0137 0.19387 0.8724
#>       H04 6.728   -0.0692 0.9249 0.08022 0.9708
#>       H05 6.697   -1.0060 1.0482 0.25099 0.8134
#>       H06 7.547   -0.1741 1.0506 0.19287 0.8812
#>       H07 8.273    0.5496 1.0509 0.17414 0.9010
#>       H08 6.832   -0.3027 0.9708 0.19796 0.8574
#>       H09 7.846   -0.2999 1.1084 0.17337 0.9109
#>       H10 7.622    1.4206 0.8438 0.07592 0.9686
#>       H11 7.208    0.1554 0.9597 0.30665 0.7100
#>       H12 7.731    1.5550 0.8404 0.04465 0.9888
#> 
#> Shukla (1972): stability variance sigma2 and s2 (heterogeneity removed), plot scale
#>  genotype  mean  sigma2       F         p       s2    F_s2      p_s2
#>       H01 6.744 3.93084 43.6780 2.154e-26  4.88526 54.2832 5.549e-27
#>       H02 7.810 0.49002  5.4449 1.376e-04  0.54579  6.0646 1.641e-04
#>       H03 7.151 0.58873  6.5418 1.811e-05  0.73893  8.2107 6.054e-06
#>       H04 6.728 0.05612  0.6235 6.821e-01  0.04266  0.4740 7.547e-01
#>       H05 6.697 1.05273 11.6976 2.275e-09  1.30700 14.5229 7.448e-10
#>       H06 7.547 0.59249  6.5835 1.678e-05  0.73033  8.1152 6.994e-06
#>       H07 8.273 0.46968  5.2189 2.098e-04  0.57669  6.4079 9.611e-05
#>       H08 6.832 0.62042  6.8939 9.517e-06  0.77482  8.6095 3.321e-06
#>       H09 7.846 0.50588  5.6211 9.916e-05  0.57070  6.3414 1.066e-04
#>       H10 7.622 0.12796  1.4218 2.205e-01  0.02763  0.3071 8.728e-01
#>       H11 7.208 1.60459 17.8296 1.775e-13  2.00071 22.2312 4.766e-14
#>       H12 7.731 0.06533  0.7260 6.051e-01 -0.05664 -0.6293 1.000e+00
#> 
#> Wricke (1962): ecovalence W (plot scale) and its share of the interaction (%)
#>  genotype  mean       W    pct
#>       H01 6.744 16.7293 36.122
#>       H02 7.810  2.3926  5.166
#>       H03 7.151  2.8039  6.054
#>       H04 6.728  0.5847  1.262
#>       H05 6.697  4.7373 10.229
#>       H06 7.547  2.8196  6.088
#>       H07 8.273  2.3079  4.983
#>       H08 6.832  2.9360  6.339
#>       H09 7.846  2.4587  5.309
#>       H10 7.622  0.8840  1.909
#>       H11 7.208  7.0366 15.193
#>       H12 7.731  0.6231  1.345
#> 
#> Kang (1993): yield-stability statistic YS
#>  genotype  mean rank adjustment adjusted_rank  sigma2 rating  ys selected
#>       H01 6.744    3         -3             0 3.93084     -8  -8    FALSE
#>       H02 7.810   10          3            13 0.49002     -8   5     TRUE
#>       H03 7.151    5         -2             3 0.58873     -8  -5    FALSE
#>       H04 6.728    2         -3            -1 0.05612      0  -1    FALSE
#>       H05 6.697    1         -3            -2 1.05273     -8 -10    FALSE
#>       H06 7.547    7          2             9 0.59249     -8   1     TRUE
#>       H07 8.273   12          3            15 0.46968     -8   7     TRUE
#>       H08 6.832    4         -3             1 0.62042     -8  -7    FALSE
#>       H09 7.846   11          3            14 0.50588     -8   6     TRUE
#>       H10 7.622    8          2            10 0.12796      0  10     TRUE
#>       H11 7.208    6         -1             5 1.60459     -8  -3    FALSE
#>       H12 7.731    9          2            11 0.06533      0  11     TRUE
#> LSD (alpha = 0.05) = 0.1978   mean YS = 0.5; selected: YS above the mean
#> 
#> Huehn (1979): rank statistics with the Z tests of Nassar and Huehn (1987)
#>  genotype  mean rank_mean    S1     Z1      p1     S2        Z2        p2
#>       H01 6.744         3 6.467 6.5572 0.01045 32.300 14.836374 0.0001173
#>       H02 7.810        10 4.400 0.1928 0.66056 14.267  0.197203 0.6569882
#>       H03 7.151         5 4.067 0.0094 0.92276 11.500  0.006199 0.9372420
#>       H04 6.728         2 2.133 3.5635 0.05906  3.067  2.796816 0.0944517
#>       H05 6.697         1 5.733 3.2685 0.07062 23.600  4.874285 0.0272597
#>       H06 7.547         7 4.533 0.3318 0.56460 13.467  0.085791 0.7695981
#>       H07 8.273        12 3.400 0.3451 0.55692  8.300  0.467083 0.4943321
#>       H08 6.832         4 3.467 0.2693 0.60377  8.267  0.475733 0.4903616
#>       H09 7.846        11 4.400 0.1928 0.66056 12.800  0.027863 0.8674315
#>       H10 7.622         8 3.000 0.9961 0.31826  6.300  1.126509 0.2885212
#>       H11 7.208         6 6.333 5.8749 0.01536 29.100 10.543679 0.0011659
#>       H12 7.731         9 2.733 1.6175 0.20345  5.367  1.532004 0.2158119
#> E(S1) = 3.972  Var(S1) = 0.9489  E(S2) = 11.92  Var(S2) = 28
#> Sum of Z1 = 23.22 (p = 0.0259)   Sum of Z2 = 36.97 (p = 0.000226); chi-square critical values: 8.21 (each Z, Bonferroni), 21.03 (sum)
st$kang
#> $table
#>    genotype     mean rank adjustment adjusted_rank     sigma2 rating  ys
#> 1       H01 6.744444    3         -3             0 3.93083687     -8  -8
#> 2       H02 7.810000   10          3            13 0.49001731     -8   5
#> 3       H03 7.151111    5         -2             3 0.58873420     -8  -5
#> 4       H04 6.727778    2         -3            -1 0.05611554      0  -1
#> 5       H05 6.697222    1         -3            -2 1.05273442     -8 -10
#> 6       H06 7.547222    7          2             9 0.59248642     -8   1
#> 7       H07 8.272778   12          3            15 0.46968154     -8   7
#> 8       H08 6.831667    4         -3             1 0.62042331     -8  -7
#> 9       H09 7.846111   11          3            14 0.50587620     -8   6
#> 10      H10 7.621667    8          2            10 0.12795665      0  10
#> 11      H11 7.208333    6         -1             5 1.60458598     -8  -3
#> 12      H12 7.731111    9          2            11 0.06533287      0  11
#>    selected
#> 1     FALSE
#> 2      TRUE
#> 3     FALSE
#> 4     FALSE
#> 5     FALSE
#> 6      TRUE
#> 7      TRUE
#> 8     FALSE
#> 9      TRUE
#> 10     TRUE
#> 11    FALSE
#> 12     TRUE
#> 
#> $lsd
#> [1] 0.1978053
#> 
#> $mean_ys
#> [1] 0.5
#> 
#> $alpha
#> [1] 0.05
#> 
```
