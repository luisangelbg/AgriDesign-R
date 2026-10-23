# List the experimental designs that can be analysed

List the experimental designs that can be analysed

## Usage

``` r
designs()
```

## Value

A data frame with the identifier of each design (the value to pass to
`design =` in
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)),
its name and the roles it needs.

## Examples

``` r
designs()
#>                 id
#> 1              crd
#> 2             rcbd
#> 3            grcbd
#> 4            latin
#> 5           graeco
#> 6           youden
#> 7              bib
#> 8              ibd
#> 9        augmented
#> 10        fact_crd
#> 11       fact_rcbd
#> 12      split_rcbd
#> 13       split_crd
#> 14      strip_rcbd
#> 15 splitsplit_rcbd
#> 16          nested
#>                                                                                 name
#> 1                                                 Completely randomised design (CRD)
#> 2                                            Randomised complete block design (RCBD)
#> 3                             Generalised RCBD (several plots per block x treatment)
#> 4                                                                       Latin square
#> 5                                                                Graeco-Latin square
#> 6                      Youden square (incomplete Latin square), intra-block analysis
#> 7  Incomplete blocks without replicate grouping (BIBD, cyclic), intra-block analysis
#> 8               Incomplete block design (lattice, alpha, BIBD), intra-block analysis
#> 9                        Augmented design (replicated checks + unreplicated entries)
#> 10                                                                  Factorial in CRD
#> 11                                                                 Factorial in RCBD
#> 12                                         Split-plot in RCBD (main plots in blocks)
#> 13                              Split-plot in CRD (main plots completely randomised)
#> 14                                                  Strip-plot (split-block) in RCBD
#> 15                                                          Split-split-plot in RCBD
#> 16                             Nested (hierarchical): second factor within the first
#>    treatments blocks row_col
#> 1           1      0   FALSE
#> 2           1      1   FALSE
#> 3           1      1   FALSE
#> 4           1      0    TRUE
#> 5           1      1    TRUE
#> 6           1      0    TRUE
#> 7           1      1   FALSE
#> 8           1      2   FALSE
#> 9           1      1   FALSE
#> 10        2-4      0   FALSE
#> 11        2-4      1   FALSE
#> 12          2      1   FALSE
#> 13          2      1   FALSE
#> 14          2      1   FALSE
#> 15          3      1   FALSE
#> 16          2      0   FALSE
```
