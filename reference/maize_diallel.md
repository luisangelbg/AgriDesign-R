# Full diallel among six maize lines

Simulated yields of a complete diallel (Griffing's method 1: parents, F1
and reciprocals) among six inbred lines, in three randomised complete
blocks. The other methods use subsets of the same data.

## Usage

``` r
maize_diallel
```

## Format

A data frame with 108 rows and 4 columns:

- Block:

  Block (1-3).

- Female:

  Female parent (P1-P6).

- Male:

  Male parent (P1-P6); equal to `Female` for the parents.

- Yield_t_ha:

  Grain yield, t/ha.

## Examples

``` r
diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block", method = 2)
#> 
#> Diallel analysis, Griffing (1956) method 2 (parents and F1)
#> Response: Yield_t_ha   6 parents, r = 3; 45 row(s) outside this method left out
#> 
#> Analysis of variance of the entries (plot scale):
#>         Df   Sum Sq Mean Sq F value    Pr(>F)    
#> Blocks   2   0.5473  0.2736   2.416 1.022e-01    
#> Entries 20 157.3162  7.8658  69.439 6.870e-25 ***
#> Error   40   4.5311  0.1133                      
#> Total   62 162.3945  2.6193                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Combining ability (scale of entry means; error = MSE / r):
#>       Df Sum Sq Mean Sq F value    Pr(>F)    
#> GCA    5  8.234 1.64675   43.61 3.716e-15 ***
#> SCA   15 44.205 2.94700   78.05 8.936e-25 ***
#> Error 40  1.510 0.03776                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Model II: GCA against SCA, F = 0.5588 on 5 and 15 df, p = 0.73
#> 
#> GCA effects:
#>  parent estimate      se       t         p
#>      P1  0.67431 0.06272 10.7518 2.291e-13
#>      P2 -0.43736 0.06272 -6.9737 2.040e-08
#>      P3  0.18681 0.06272  2.9786 4.902e-03
#>      P4 -0.55653 0.06272 -8.8739 5.385e-11
#>      P5  0.18847 0.06272  3.0052 4.566e-03
#>      P6 -0.05569 0.06272 -0.8881 3.798e-01
#> 
#> SCA effects:
#>  parent1 parent2 estimate     se       t         p
#>       P1      P1  -2.3704 0.1422 -16.666 1.334e-19
#>       P1      P2   0.4880 0.1722   2.833 7.188e-03
#>       P1      P3   1.0205 0.1722   5.925 6.033e-07
#>       P1      P4   0.9938 0.1722   5.770 9.957e-07
#>       P1      P5   1.2988 0.1722   7.541 3.351e-09
#>       P1      P6   0.9396 0.1722   5.455 2.749e-06
#>       P2      P2  -2.0770 0.1422 -14.604 1.268e-17
#>       P2      P3   1.3788 0.1722   8.005 7.782e-10
#>       P2      P4   0.5655 0.1722   3.283 2.138e-03
#>       P2      P5   0.7638 0.1722   4.434 7.031e-05
#>       P2      P6   0.9580 0.1722   5.562 1.950e-06
#>       P3      P3  -2.3787 0.1422 -16.725 1.179e-19
#>       P3      P4   1.5446 0.1722   8.968 4.055e-11
#>       P3      P5  -0.1904 0.1722  -1.105 2.757e-01
#>       P3      P6   1.0038 0.1722   5.828 8.252e-07
#>       P4      P4  -2.2420 0.1422 -15.764 9.273e-19
#>       P4      P5   0.5763 0.1722   3.346 1.793e-03
#>       P4      P6   0.8038 0.1722   4.667 3.402e-05
#>       P5      P5  -1.8120 0.1422 -12.741 1.158e-15
#>       P5      P6   1.1755 0.1722   6.825 3.294e-08
#>       P6      P6  -2.4404 0.1422 -17.158 4.787e-20
#> 
#> Standard errors of differences:
#>   difference      se
#>    g_i - g_j 0.09716
#>  s_ii - s_jj 0.19432
#>  s_ij - s_ik 0.25706
#>  s_ij - s_kl 0.23799
#> 
#> Model I (fixed effects):
#>          quantity estimate
#>  sum(g^2)/(p - 1)   0.2011
#>       sum(s^2)/df   2.9092
#> 
#> Model II (random effects), entry-mean scale:
#>                  component estimate       se
#>                 sigma2 GCA -0.16253 0.187195
#>                 sigma2 SCA  2.90924 1.076125
#>  sigma2 error (entry mean)  0.03776 0.008443
#> 
#> Genetic parameters (inbreeding of the parents F = 1):
#>                                               parameter estimate
#>                                                sigma2_A  -0.3251
#>                                                sigma2_D   2.9092
#>             dominance ratio sqrt(2 sigma2_D / sigma2_A)         
#>  Baker ratio 2 sigma2_GCA / (2 sigma2_GCA + sigma2_SCA)         
```
