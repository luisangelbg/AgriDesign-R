# AUDPC and AUDPS of every plot of a trial

Takes the assessments in long format (one row per plot and date), sorts
them by time within each plot and returns one row per plot with its
areas, ready for
[`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md).

## Usage

``` r
audpc_by_unit(data, response, time, units, ymax = 100)
```

## Arguments

- data:

  A data frame in long format.

- response:

  Column with the disease intensity.

- time:

  Column with the time of assessment (numbers or
  [Date](https://rdrr.io/r/base/Dates.html)).

- units:

  Column(s) that identify a plot, for instance
  `c("Block", "Treatment")`; they are kept in the result.

- ymax:

  Largest possible intensity (100 for percentages, 1 for proportions,
  the top of a rating scale).

## Value

A data frame with the unit columns, the number of assessments (`n`), the
time spanned (`duration`) and the areas `audpc`, `audps`, their relative
(`audpc_rel`, `audps_rel`) and standardised (`audpc_std`, `audps_std`)
forms.

## Examples

``` r
a <- audpc_by_unit(blight_progress, "Severity_pct", "Day", c("Block", "Treatment"))
head(a)
#>   Block  Treatment n duration   audpc  audps audpc_rel audps_rel audpc_std
#> 1     1  Untreated 7       42 3484.25 3917.9 0.8295833 0.7995714  82.95833
#> 2     1    Contact 7       42 2480.80 2856.7 0.5906667 0.5830000  59.06667
#> 3     1   Systemic 7       42 1948.10 2294.6 0.4638333 0.4682857  46.38333
#> 4     1 Alternated 7       42 1693.65 2005.5 0.4032500 0.4092857  40.32500
#> 5     2  Untreated 7       42 3259.20 3659.6 0.7760000 0.7468571  77.60000
#> 6     2    Contact 7       42 2489.55 2849.7 0.5927500 0.5815714  59.27500
#>   audps_std
#> 1  79.95714
#> 2  58.30000
#> 3  46.82857
#> 4  40.92857
#> 5  74.68571
#> 6  58.15714
fit <- design_anova(a, "audps", "Treatment", block = "Block")
compare_means(fit, "Treatment")
#> 
#> Tukey HSD (Tukey-Kramer), alpha = 0.05
#> Response: audps   compared: Treatment
#> Error: Residual error (MS = 20409, df = 9)
#> Minimum significant difference = 315.4
#> 
#>   Treatment n mean    se raw_mean     sd group
#>   Untreated 4 3591 71.43     3591 252.38     a
#>     Contact 4 2898 71.43     2898  87.68     b
#>    Systemic 4 2105 71.43     2105 159.93     c
#>  Alternated 4 1788 71.43     1788 225.57     d
#> 
#> Means sharing a letter do not differ.
```
