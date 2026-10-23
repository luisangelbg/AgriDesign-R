# Rank-based and unequal-variance alternatives to the ANOVA

When the assumptions fail, compares the treatments without assuming
normality (rank-based tests) or without assuming equal variances
(Welch). The method follows the design of `fit`:

## Usage

``` r
nonparametric(
  fit,
  method = NULL,
  posthoc = c("conover", "nemenyi"),
  p_adjust = "holm",
  alpha = 0.05
)
```

## Arguments

- fit:

  An `"ad_anova"` object from
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

- method:

  One of the methods above; `NULL` picks the first one that fits the
  design (Kruskal-Wallis, Friedman, Scheirer-Ray-Hare or ART).

- posthoc:

  For Friedman, `"conover"` (default) or `"nemenyi"`.

- p_adjust:

  Adjustment of the pairwise p-values, as in
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- alpha:

  Significance level of the checks.

## Value

An object of class `"ad_nonparametric"` with the test (`test`), the
table of terms for SRH and ART (`table`), the pairwise comparisons
(`pairs`) and the groups with their letters (`groups`).

## Details

- `"kruskal"`: Kruskal-Wallis with Dunn's pairwise tests (one factor, no
  blocks); `"kruskal_mw"` uses pairwise Mann-Whitney tests instead.

- `"welch"`: Welch's one-way ANOVA with Games-Howell pairwise tests (one
  factor, no blocks; normal data with unequal variances).

- `"friedman"`: Friedman's test with Conover's or Nemenyi's pairwise
  tests (one factor in complete blocks, one plot per block and
  treatment).

- `"srh"`: Scheirer-Ray-Hare extension of Kruskal-Wallis (two factors,
  no blocks).

- `"art"`: aligned rank transform ANOVA (any factorial or blocked design
  tested against the residual error), with Tukey comparisons on the
  aligned ranks of the first treatment factor.

Pairwise p-values of Dunn, Mann-Whitney and Conover are adjusted with
`p_adjust`; Games-Howell, Nemenyi and Tukey already control the
family-wise error. Levels that share a letter do not differ.

## Examples

``` r
fit <- design_anova(bean_varieties, "Yield_kg_plot", "Variety")
nonparametric(fit)
#> 
#> Kruskal-Wallis rank sum test
#> Response: Yield_kg_plot   factor: Variety
#> H = 17.63, df = 4, p = 0.001458
#> epsilon2 = 0.9086
#> 
#> Pairwise comparisons: Dunn (p adjusted by holm)
#>         level_i        level_j statistic         p    p_adj different
#>    Negro Jamapa   Flor de Mayo   -1.9721 0.0485951 0.340166     FALSE
#>    Negro Jamapa Pinto Saltillo   -1.0757 0.2820589 1.000000     FALSE
#>    Negro Jamapa           Bayo    0.8367 0.4027837 1.000000     FALSE
#>    Negro Jamapa       Azufrado   -2.8685 0.0041236 0.037112      TRUE
#>    Flor de Mayo Pinto Saltillo    0.8964 0.3700277 1.000000     FALSE
#>    Flor de Mayo           Bayo    2.8088 0.0049728 0.039783      TRUE
#>    Flor de Mayo       Azufrado   -0.8964 0.3700277 1.000000     FALSE
#>  Pinto Saltillo           Bayo    1.9124 0.0558293 0.340166     FALSE
#>  Pinto Saltillo       Azufrado   -1.7928 0.0729980 0.364990     FALSE
#>            Bayo       Azufrado   -3.7052 0.0002112 0.002112      TRUE
#> 
#> Groups (levels sharing a letter do not differ at alpha = 0.05):
#>           level n median  mean mean_rank group
#>        Azufrado 4  2.820 2.812     18.25     a
#>    Flor de Mayo 4  2.595 2.598     14.50    ab
#>  Pinto Saltillo 4  2.425 2.420     10.75   abc
#>    Negro Jamapa 4  2.115 2.132      6.25    bc
#>            Bayo 4  1.920 1.917      2.75     c

fit <- design_anova(disease_scores, "Severity_score", "Fungicide", block = "Block")
nonparametric(fit, posthoc = "nemenyi")
#> 
#> Friedman rank sum test
#> Response: Severity_score   factor: Fungicide
#> chisq = 14.04, df = 3, p = 0.002851
#> Kendall_W = 0.936
#> Iman-Davenport F = 58.5 on 3 and 12 df, p = 1.959e-07
#> 
#> Pairwise comparisons: Nemenyi
#>       level_i      level_j statistic        p    p_adj different
#>       Control     Mancozeb     3.118 0.121943 0.121943     FALSE
#>       Control Azoxystrobin     5.196 0.001363 0.001363      TRUE
#>       Control       Copper     2.078 0.455933 0.455933     FALSE
#>      Mancozeb Azoxystrobin     2.078 0.455933 0.455933     FALSE
#>      Mancozeb       Copper     1.039 0.883087 0.883087     FALSE
#>  Azoxystrobin       Copper     3.118 0.121943 0.121943     FALSE
#> 
#> Groups (levels sharing a letter do not differ at alpha = 0.05):
#>         level n median mean mean_rank group
#>       Control 5      8  7.8       4.0     a
#>        Copper 5      6  5.8       2.8    ab
#>      Mancozeb 5      5  4.8       2.2    ab
#>  Azoxystrobin 5      3  2.6       1.0     b
```
