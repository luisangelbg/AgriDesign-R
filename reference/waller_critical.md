# Critical t of the Waller-Duncan k-ratio test

The Bayesian critical value t_B(K, F, q, f) of Waller and Duncan (1969),
computed by numerical integration of the posterior expected losses (see
[`waller_duncan()`](https://luisangelbg.github.io/AgriDesign-R/reference/waller_duncan.md)).
Two means differ when their difference exceeds `t_B * SE`, the Bayes
least significant difference. The value decreases as the observed F
grows: with a small F the means are shrunk towards each other and few
differences are declared.

## Usage

``` r
waller_critical(K = 100, q, f, F)
```

## Arguments

- K:

  Ratio of the seriousness of a type I error to that of a type II error
  (50, 100 and 500 roughly match alpha = 0.10, 0.05 and 0.01).

- q, f:

  Degrees of freedom of the treatments and of the error.

- F:

  Observed F ratio of the treatments.

## Value

A numeric vector of critical values (arguments are recycled).

## References

Waller, R. A. and Duncan, D. B. (1969). A Bayes rule for the symmetric
multiple comparisons problem. *Journal of the American Statistical
Association* 64, 1484-1503.

Waller, R. A. and Kemp, K. E. (1975). Computations of Bayesian t-values
for multiple comparisons. *Journal of Statistical Computation and
Simulation* 4, 169-171.

## Examples

``` r
waller_critical(100, q = 5, f = 24, F = 14.37)
#> [1] 1.918732
waller_critical(100, q = 4, f = 20, F = c(2, 5, 10, 100))
#> [1] 2.508500 2.143551 1.991624 1.873337
```
