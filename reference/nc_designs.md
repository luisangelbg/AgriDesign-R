# North Carolina mating designs

Simulated data of the three North Carolina designs, each in two sets
with two replicates per set.

## Usage

``` r
nc_designs
```

## Format

A list of three data frames:

- nc1:

  Design I: `Set`, `Rep`, `Male` (three per set), `Female` (four per
  male, labelled within the male) and `Yield_t_ha`; 48 rows.

- nc2:

  Design II: `Set`, `Rep`, `Male` and `Female` (four of each per set)
  and `Yield_t_ha`; 64 rows.

- nc3:

  Design III: `Set`, `Rep`, `Male` (six F2 plants per set), `Tester`
  (`P1`, `P2`) and `Yield_t_ha`; 48 rows.

## Examples

``` r
nc_design(nc_designs$nc1, "Yield_t_ha", "Male", "Female", block = "Rep", set = "Set", design = 1)
#> 
#> North Carolina design I (females nested within males)
#> Response: Yield_t_ha   sets = 2, blocks per set = 2, males per set = 3, females per male = 4
#> 
#>                Df Sum Sq Mean Sq F value    Pr(>F)    
#> Sets            1  2.881 2.88120                      
#> Blocks(sets)    2  2.812 1.40622  16.026 5.079e-05 ***
#> Males(sets)     4  9.413 2.35330   2.764 5.945e-02   .
#> Females(males) 18 15.325 0.85140   9.703 1.096e-06 ***
#> Error          22  1.930 0.08775                      
#> Total          47 32.363 0.68857                      
#> 
#> Tested against another mean square:
#>   Males(sets) vs Females(males)
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expected mean squares:
#>   Males(sets)             Ve + 2 Vf + 8 Vm
#>   Females(males)          Ve + 2 Vf
#>   Error                   Ve
#> 
#> Variance components:
#>                  component estimate     se
#>                 Vm (males)   0.1877 0.2110
#>  Vf (females within males)   0.3818 0.1425
#> 
#> Genetic parameters (inbreeding of the parents F = 0):
#>                                    parameter estimate    se
#>                                  Cov HS = Vm   0.1877 0.211
#>                             Cov FS = Vm + Vf   0.5696 0.234
#>                                     sigma2_A   0.7510 0.844
#>                                     sigma2_D   0.7764 1.095
#>  dominance ratio sqrt(2 sigma2_D / sigma2_A)   1.4379      
```
