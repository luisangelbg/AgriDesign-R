# Diallel analysis by Griffing's methods

Griffing's (1956) analysis of a diallel cross by any of the four
methods: 1, parents, F1 and reciprocals; 2, parents and one F1 per pair;
3, F1 and reciprocals; 4, one F1 per pair. The analysis of variance of
general (GCA) and specific (SCA) combining ability, and of reciprocal
effects (split into maternal and non-maternal, Cockerham 1963) in
methods 1 and 3; the effects with their exact standard errors (Model I),
the standard errors of differences, and the variance components of Model
II with the additive and dominance variances.

## Usage

``` r
diallel(
  data,
  response,
  parent1,
  parent2,
  block = NULL,
  method = 1,
  inbreeding = 1
)
```

## Arguments

- data:

  A data frame in long format.

- response:

  Name of the numeric response column.

- parent1, parent2:

  Names of the female and male parent columns (equal for the parents
  themselves).

- block:

  Name of the block column (randomised complete blocks); when `NULL` the
  entries are taken as completely randomised.

- method:

  Griffing's method, 1 to 4.

- inbreeding:

  Inbreeding coefficient of the parents, used to turn the Model II
  components into additive and dominance variances (1 for inbred lines).

## Value

An object of class `"ad_diallel"`.

## Details

The combining-ability analysis is on the scale of entry means: the error
is the error mean square of the trial divided by the number of
replicates. Rows that do not belong to the method are left out (parents
in methods 3 and 4; reciprocals in methods 2 and 4, where the cross with
the first parent as female is kept), so one data set serves all methods.
The effects are the least-squares estimates under Griffing's
restrictions (sum of GCA effects zero; SCA effects summing to zero for
every parent, the selfs counted twice in method 2), which reproduce his
closed formulas.

## References

Griffing, B. (1956) Australian Journal of Biological Sciences 9,
463-493. Cockerham, C. C. (1963) In Statistical Genetics and Plant
Breeding, NAS-NRC 982, 53-94.

## Examples

``` r
dl <- diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block", method = 1)
dl
#> 
#> Diallel analysis, Griffing (1956) method 1 (parents, F1 and reciprocals)
#> Response: Yield_t_ha   6 parents, r = 3
#> 
#> Analysis of variance of the entries (plot scale):
#>          Df  Sum Sq Mean Sq F value    Pr(>F)    
#> Blocks    2   1.124  0.5619   4.919 1.003e-02   *
#> Entries  35 188.652  5.3901  47.187 1.359e-36 ***
#> Error    70   7.996  0.1142                      
#> Total   107 197.772  1.8483                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Combining ability (scale of entry means; error = MSE / r):
#>                Df  Sum Sq Mean Sq F value    Pr(>F)    
#> GCA             5 10.9555 2.19110 57.5453 1.914e-23 ***
#> SCA            15 51.1168 3.40779 89.4993 1.730e-39 ***
#> Reciprocal     15  0.8117 0.05411  1.4212 1.619e-01    
#>   Maternal      5  0.5603 0.11206  2.9429 1.810e-02   *
#>   Non-maternal 10  0.2514 0.02514  0.6604 7.568e-01    
#> Error          70  2.6653 0.03808                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Model II: GCA against SCA, F = 0.643 on 5 and 15 df, p = 0.671 (approximate in method 1)
#> 
#> GCA effects:
#>  parent estimate      se       t         p
#>      P1  0.63620 0.05142  12.372 2.752e-19
#>      P2 -0.39435 0.05142  -7.669 7.429e-11
#>      P3  0.14315 0.05142   2.784 6.903e-03
#>      P4 -0.53296 0.05142 -10.365 8.829e-16
#>      P5  0.21037 0.05142   4.091 1.132e-04
#>      P6 -0.06241 0.05142  -1.214 2.290e-01
#> 
#> SCA effects:
#>  parent1 parent2 estimate     se        t         p
#>       P1      P1 -2.61648 0.1626 -16.0906 3.229e-25
#>       P1      P2  0.15907 0.1173   1.3566 1.793e-01
#>       P1      P3  0.67657 0.1173   5.7699 1.986e-07
#>       P1      P4  0.48602 0.1173   4.1448 9.388e-05
#>       P1      P5  0.88269 0.1173   7.5276 1.351e-10
#>       P1      P6  0.41213 0.1173   3.5147 7.765e-04
#>       P2      P2 -2.48537 0.1626 -15.2843 5.374e-24
#>       P2      P3  0.97713 0.1173   8.3331 4.448e-12
#>       P2      P4  0.07991 0.1173   0.6815 4.978e-01
#>       P2      P5  0.61991 0.1173   5.2866 1.348e-06
#>       P2      P6  0.64935 0.1173   5.5377 5.018e-07
#>       P3      P3 -2.61370 0.1626 -16.0735 3.424e-25
#>       P3      P4  1.06407 0.1173   9.0745 1.930e-13
#>       P3      P5 -0.53259 0.1173  -4.5420 2.263e-05
#>       P3      P6  0.42852 0.1173   3.6545 4.942e-04
#>       P4      P4 -2.61148 0.1626 -16.0599 3.589e-25
#>       P4      P5  0.46519 0.1173   3.9672 1.736e-04
#>       P4      P6  0.51630 0.1173   4.4030 3.750e-05
#>       P5      P5 -2.17815 0.1626 -13.3950 5.402e-21
#>       P5      P6  0.74296 0.1173   6.3361 1.980e-08
#>       P6      P6 -2.74926 0.1626 -16.9072 2.032e-26
#> 
#> Reciprocal effects r_ij = (x_ij - x_ji) / 2:
#>  parent1 parent2  estimate    se        t       p
#>       P1      P2  0.001667 0.138  0.01208 0.99040
#>       P1      P3  0.103333 0.138  0.74891 0.45642
#>       P1      P4  0.200000 0.138  1.44950 0.15166
#>       P1      P5  0.110000 0.138  0.79723 0.42802
#>       P1      P6  0.250000 0.138  1.81188 0.07429
#>       P2      P3  0.080000 0.138  0.57980 0.56391
#>       P2      P4  0.096667 0.138  0.70059 0.48588
#>       P2      P5 -0.243333 0.138 -1.76356 0.08217
#>       P2      P6 -0.050000 0.138 -0.36238 0.71816
#>       P3      P4  0.178333 0.138  1.29247 0.20044
#>       P3      P5  0.041667 0.138  0.30198 0.76356
#>       P3      P6  0.303333 0.138  2.19841 0.03123
#>       P4      P5 -0.256667 0.138 -1.86019 0.06706
#>       P4      P6 -0.051667 0.138 -0.37445 0.70920
#>       P5      P6  0.095000 0.138  0.68851 0.49341
#> 
#> Maternal effects:
#>  parent estimate      se       t       p
#>      P1  0.11083 0.05142  2.1554 0.03457
#>      P2 -0.01972 0.05142 -0.3835 0.70248
#>      P3  0.05667 0.05142  1.1020 0.27424
#>      P4 -0.13056 0.05142 -2.5389 0.01334
#>      P5  0.07389 0.05142  1.4369 0.15519
#>      P6 -0.09111 0.05142 -1.7718 0.08077
#> 
#> Standard errors of differences:
#>   difference      se
#>    g_i - g_j 0.07966
#>  s_ii - s_jj 0.22532
#>  s_ij - s_ik 0.17813
#>  s_ij - s_kl 0.15932
#>  r_ij - r_kl 0.19513
#> 
#> Model I (fixed effects):
#>          quantity estimate
#>  sum(g^2)/(p - 1) 0.179419
#>       sum(s^2)/df 1.684855
#>       sum(r^2)/df 0.008019
#> 
#> Model II (random effects), entry-mean scale:
#>                  component  estimate       se
#>                 sigma2 GCA -0.092332 0.152991
#>                 sigma2 SCA  1.956606 0.722534
#>          sigma2 reciprocal  0.008019 0.010391
#>            sigma2 maternal  0.007243 0.005980
#>        sigma2 non-maternal -0.006466 0.006478
#>  sigma2 error (entry mean)  0.038076 0.006436
#> 
#> Genetic parameters (inbreeding of the parents F = 1):
#>                                               parameter estimate
#>                                                sigma2_A  -0.1847
#>                                                sigma2_D   1.9566
#>             dominance ratio sqrt(2 sigma2_D / sigma2_A)         
#>  Baker ratio 2 sigma2_GCA / (2 sigma2_GCA + sigma2_SCA)         
diallel(maize_diallel, "Yield_t_ha", "Female", "Male", block = "Block", method = 4)$gca
#>   parent    estimate         se          t            p
#> 1     P1  0.74944444 0.09599325  7.8072617 1.667550e-08
#> 2     P2 -0.50888889 0.09599325 -5.3012985 1.217514e-05
#> 3     P3  0.26611111 0.09599325  2.7721856 9.788312e-03
#> 4     P4 -0.54555556 0.09599325 -5.6832698 4.305838e-06
#> 5     P5 -0.01555556 0.09599325 -0.1620484 8.724315e-01
#> 6     P6  0.05444444 0.09599325  0.5671695 5.751184e-01
```
