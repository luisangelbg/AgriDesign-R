# List the designs that can be analysed or randomised

`designs()` lists the 16 designs that
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
analyses; `designs("layout")` lists the 16 designs that
[`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)
randomises, with the arguments each one needs and the design to give to
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
once the responses are in the field book.

## Usage

``` r
designs(type = c("analysis", "layout"))
```

## Arguments

- type:

  `"analysis"` (designs of
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md))
  or `"layout"` (designs of
  [`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)).

## Value

A data frame. For `"analysis"`: the identifier of each design (the value
to pass to `design =` in
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)),
its name, the number of treatment factors and block columns, and whether
it needs rows and columns. For `"layout"`: the identifier (the first
argument of
[`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)),
its name, the arguments it needs (in parentheses when optional), its
conditions and the design of
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
that analyses its field book.

## Details

Designs analysed by
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
(`design =`):

|  |  |  |
|----|----|----|
| id | Design | Roles |
| `crd` | Completely randomised | treatment |
| `rcbd` | Randomised complete blocks | treatment, block |
| `grcbd` | Generalised RCBD (replicates within blocks) | treatment, block |
| `latin` | Latin square | treatment, row, col |
| `graeco` | Graeco-Latin square | treatment, row, col, block (Greek letters) |
| `youden` | Youden square | treatment, row, col |
| `bib` | Incomplete blocks (BIBD, cyclic) | treatment, block |
| `ibd` | Resolvable incomplete blocks (lattice, alpha) | treatment, block = c(rep, block) |
| `augmented` | Augmented (checks and new entries) | treatment, block |
| `fact_crd` | Factorial in CRD (2 to 4 factors) | treatment = c(A, B, ...) |
| `fact_rcbd` | Factorial in RCBD (2 to 4 factors) | treatment = c(A, B, ...), block |
| `split_rcbd` | Split-plot in RCBD | treatment = c(main, sub), block |
| `split_crd` | Split-plot in CRD | treatment = c(main, sub), block (replicate) |
| `strip_rcbd` | Strip-plot (split-block) | treatment = c(A, B), block |
| `splitsplit_rcbd` | Split-split-plot in RCBD | treatment = c(main, sub, sub-sub), block |
| `nested` | Nested (hierarchical) | treatment = c(outer, inner) |

Designs randomised by
[`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)
(first argument; `C` is optional in factorials, and the conditions on
the number of treatments and the block size are in the `conditions`
column of `designs("layout")`):

|  |  |  |  |
|----|----|----|----|
| id | Design | Arguments | Analyse as |
| `crd` | Completely randomised | `trt`, `r` | `crd` |
| `rcbd` | Randomised complete blocks | `trt`, `r` | `rcbd` |
| `latin` | Latin square | `trt` | `latin` |
| `graeco` | Graeco-Latin square | `trt`, `greek` | `graeco` |
| `youden` | Youden square | `trt`, `k` | `youden` |
| `bib` | Balanced incomplete blocks | `trt`, `k` | `bib` |
| `cyclic` | Cyclic incomplete blocks | `trt`, `k`, `r` | `bib` |
| `ibd` | Resolvable incomplete blocks | `trt`, `k`, `r` | `ibd` |
| `alpha` | Alpha design | `trt`, `k`, `r` | `ibd` |
| `lattice` | Square lattice | `trt`, `k`, `r` | `ibd` |
| `augmented` | Augmented | `checks`, `entries`, `r` | `augmented` |
| `fact_crd` | Factorial in CRD | `A`, `B`, `C`, `r` | `fact_crd` |
| `fact_rcbd` | Factorial in RCBD | `A`, `B`, `C`, `r` | `fact_rcbd` |
| `split_rcbd` | Split-plot in RCBD | `A`, `B`, `r` | `split_rcbd` |
| `strip` | Strip-plot | `A`, `B`, `r` | `strip_rcbd` |
| `splitsplit` | Split-split-plot | `A`, `B`, `C`, `r` | `splitsplit_rcbd` |

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
designs("layout")
#>            id
#> 1         crd
#> 2        rcbd
#> 3       latin
#> 4    fact_crd
#> 5   fact_rcbd
#> 6  split_rcbd
#> 7       strip
#> 8  splitsplit
#> 9         ibd
#> 10  augmented
#> 11     graeco
#> 12     youden
#> 13        bib
#> 14     cyclic
#> 15      alpha
#> 16    lattice
#>                                                                        name
#> 1                                        Completely randomised design (CRD)
#> 2                                   Randomised complete block design (RCBD)
#> 3                                                              Latin square
#> 4                                                          Factorial in CRD
#> 5                                                         Factorial in RCBD
#> 6                                                        Split-plot in RCBD
#> 7                                          Strip-plot (split-block) in RCBD
#> 8                                                  Split-split-plot in RCBD
#> 9  Resolvable incomplete blocks (random blocks of size k within replicates)
#> 10          Augmented design (checks in every block + unreplicated entries)
#> 11                                                      Graeco-Latin square
#> 12                                                            Youden square
#> 13                                  Balanced incomplete block design (BIBD)
#> 14                                           Cyclic incomplete block design
#> 15                              Alpha design (resolvable incomplete blocks)
#> 16                                                           Square lattice
#>                 needs
#> 1              trt, r
#> 2              trt, r
#> 3                 trt
#> 4         A, B (C), r
#> 5         A, B (C), r
#> 6             A, B, r
#> 7             A, B, r
#> 8          A, B, C, r
#> 9           trt, k, r
#> 10 checks, entries, r
#> 11        trt (greek)
#> 12             trt, k
#> 13             trt, k
#> 14          trt, k, r
#> 15          trt, k, r
#> 16          trt, k, r
#>                                                                     conditions
#> 1                                                                             
#> 2                                                                             
#> 3                                        rows = columns = number of treatments
#> 4                                               all combinations of the levels
#> 5                                all combinations of the levels in every block
#> 6                                              A on main plots, B on sub-plots
#> 7                                 A in horizontal strips, B in vertical strips
#> 8                          A on main plots, B on sub-plots, C on sub-sub-plots
#> 9  treatments a multiple of k; each replicate split at random into blocks of k
#> 10                                      checks in every block, each entry once
#> 11                                         odd number of treatments, or 4 or 8
#> 12                                       k rows, as many columns as treatments
#> 13                             every pair of treatments together equally often
#> 14                                                           r a multiple of k
#> 15                          treatments a multiple of k and k <= treatments / k
#> 16                                                            treatments = k^2
#>         analyse_as
#> 1              crd
#> 2             rcbd
#> 3            latin
#> 4         fact_crd
#> 5        fact_rcbd
#> 6       split_rcbd
#> 7       strip_rcbd
#> 8  splitsplit_rcbd
#> 9              ibd
#> 10       augmented
#> 11          graeco
#> 12          youden
#> 13             bib
#> 14             bib
#> 15             ibd
#> 16             ibd

# the layout reports the call that analyses it
make_design("strip", A = c("Drip", "Furrow"), B = c("V1", "V2", "V3"), r = 3, seed = 1)$analysis
#> [1] "design_anova(data, \"response\", c(\"A\", \"B\"), block = \"Block\", design = \"strip_rcbd\")"
```
