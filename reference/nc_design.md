# North Carolina designs I, II and III

Analysis of the mating designs of Comstock and Robinson (1948, 1952):

- Design I:

  each male crossed with its own group of females (females nested within
  males).

- Design II:

  every male crossed with every female of its set.

- Design III:

  F2 (or other) plants used as males, each crossed with the two inbred
  parents of the population (the testers, given in the `female` column).

Males (and females) are nested within sets, and blocks within sets. When
a plot holds several plants, the plot error (blocks x crosses within
sets) is separated from the variation within plots. The table gives the
expected mean squares; every term is tested against the mean square that
differs from it only by its own component, and the components are their
moment estimators with standard errors from Var(MS) = 2 MS^2 / df.

## Usage

``` r
nc_design(
  data,
  response,
  male,
  female,
  block,
  set = NULL,
  design = 1,
  inbreeding = 0
)
```

## Arguments

- data:

  A data frame in long format.

- response:

  Name of the numeric response column.

- male, female:

  Names of the male and female columns (design III: the female column
  holds the two testers).

- block:

  Name of the block (replicate) column, nested within sets.

- set:

  Optional name of the set column.

- design:

  1, 2 or 3.

- inbreeding:

  Inbreeding coefficient of the parents (designs I and II).

## Value

An object of class `"ad_nc_design"`: the analysis of variance with
expected mean squares, the variance components and the genetic
parameters.

## Details

Genetic interpretation without epistasis: design I,
`sigma2_A = 4 Vm / (1 + F)` and `sigma2_D = 4 (Vf - Vm) / (1 + F)^2`;
design II, `sigma2_A = 2 (Vm + Vf) / (1 + F)` and
`sigma2_D = 4 Vmf / (1 + F)^2`; design III (F2 reference population,
gene frequency 1/2), `sigma2_A = 4 Vm = 2 (MSm - MSe) / r` and
`sigma2_D = Vmt = (MSmt - MSe) / r`, because the sum and the difference
of the two backcross families of a plant have variances sigma2_A and 2
sigma2_D. The average degree of dominance is
`sqrt(2 sigma2_D / sigma2_A)`.

## References

Comstock, R. E. and Robinson, H. F. (1948) Biometrics 4, 254-266.
Comstock, R. E. and Robinson, H. F. (1952) In Heterosis, Iowa State
College Press, 494-516.

## Examples

``` r
nc_design(nc_designs$nc2, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 2)
#> 
#> North Carolina design II (males x females)
#> Response: Yield_t_ha   sets = 2, blocks per set = 2, males per set = 4, females per set = 4
#> 
#>                       Df Sum Sq Mean Sq F value    Pr(>F)    
#> Sets                   1  9.602  9.6023                      
#> Blocks(sets)           2  2.489  1.2447  11.938 0.0001534 ***
#> Males(sets)            6 15.595  2.5992   8.520 0.0001777 ***
#> Females(sets)          6  9.972  1.6620   5.448 0.0022916  **
#> Males x Females(sets) 18  5.491  0.3051   2.926 0.0045367  **
#> Error                 30  3.128  0.1043                      
#> Total                 63 46.278  0.7346                      
#> 
#> Tested against another mean square:
#>   Males(sets) vs Males x Females(sets)
#>   Females(sets) vs Males x Females(sets)
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expected mean squares:
#>   Males(sets)             Ve + 2 Vmf + 8 Vm
#>   Females(sets)           Ve + 2 Vmf + 8 Vf
#>   Males x Females(sets)   Ve + 2 Vmf
#>   Error                   Ve
#> 
#> Variance components:
#>              component estimate      se
#>             Vm (males)   0.2868 0.18801
#>           Vf (females)   0.1696 0.12061
#>  Vmf (males x females)   0.1004 0.05259
#> 
#> Genetic parameters (inbreeding of the parents F = 0):
#>                                    parameter estimate     se
#>                             sigma2_A (males)   1.1471 0.7521
#>                           sigma2_A (females)   0.6785 0.4825
#>                            sigma2_A (pooled)   0.9128 0.4482
#>                                     sigma2_D   0.4016 0.2104
#>  dominance ratio sqrt(2 sigma2_D / sigma2_A)   0.9380       
nc_design(nc_designs$nc3, "Yield_t_ha", "Male", "Tester", block = "Rep", set = "Set", design = 3)
#> 
#> North Carolina design III (backcrosses to both parents)
#> Response: Yield_t_ha   sets = 2, blocks per set = 2, males per set = 6, testers = 2
#> 
#>                       Df    Sum Sq   Mean Sq  F value    Pr(>F)    
#> Sets                   1  0.003502  0.003502                       
#> Blocks(sets)           2  2.308321  1.154160  23.8233 3.124e-06 ***
#> Testers                1 52.689252 52.689252 524.2864 5.697e-10 ***
#> Testers x Sets         1  0.019602  0.019602   0.1951 6.681e-01    
#> Males(sets)           10  2.930721  0.293072   6.0494 2.168e-04 ***
#> Males x Testers(sets) 10  1.004971  0.100497   2.0744 7.393e-02   .
#> Error                 22  1.065829  0.048447                       
#> Total                 47 60.022198  1.277068                       
#> 
#> Tested against another mean square:
#>   Testers vs Males x Testers(sets)
#>   Testers x Sets vs Males x Testers(sets)
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expected mean squares:
#>   Testers                 Ve + 2 Vmt + fixed testers
#>   Males(sets)             Ve + 4 Vm
#>   Males x Testers(sets)   Ve + 2 Vmt
#>   Error                   Ve
#> 
#> Variance components:
#>              component estimate      se
#>             Vm (males)  0.06116 0.03297
#>  Vmt (males x testers)  0.02603 0.02363
#> 
#> Genetic parameters (F2 reference population):
#>                                    parameter estimate      se
#>                                     sigma2_A  0.24463 0.13188
#>                                     sigma2_D  0.02603 0.02363
#>  dominance ratio sqrt(2 sigma2_D / sigma2_A)  0.46128        
```
