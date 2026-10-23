# Area under the disease progress curve and stairs

`audpc()` integrates the disease progress curve with trapezoids (Shaner
and Finney 1977); `audps()` gives the area under the disease progress
stairs of Simko and Piepho (2012), which weighs the first and the last
assessment as much as the others:
`AUDPS = AUDPC + (y1 + yn) / 2 * D / (n - 1)`, with `D` the time from
the first to the last assessment.

## Usage

``` r
audpc(
  evaluation,
  dates,
  type = c("absolute", "relative", "standardized"),
  ymax = 100
)

audps(
  evaluation,
  dates,
  type = c("absolute", "relative", "standardized"),
  ymax = 100
)
```

## Arguments

- evaluation:

  Disease intensities: a vector (one unit) or a matrix or data frame
  with one row per unit and one column per assessment.

- dates:

  Times of the assessments (days or
  [Date](https://rdrr.io/r/base/Dates.html) values), strictly
  increasing.

- type:

  `"absolute"`, `"relative"` or `"standardized"`.

- ymax:

  Largest possible intensity (100 for percentages, 1 for proportions,
  the top of a rating scale).

## Value

A numeric vector with one area per unit.

## Details

`type = "relative"` divides the area by the largest possible one, so it
lies between 0 and 1: `ymax * D` for the AUDPC and
`ymax * D * n / (n - 1)` for the AUDPS, whose stairs are `n / (n - 1)`
times as wide as the curve. `type = "standardized"` divides by the width
only and gives the mean intensity over the epidemic, in the units of the
assessments. Missing assessments of a unit are skipped (its area then
spans its own assessments only).

## References

Shaner, G. and Finney, R. E. (1977). The effect of nitrogen
fertilization on the expression of slow-mildewing resistance in Knox
wheat. *Phytopathology* 67, 1051-1056.

Simko, I. and Piepho, H.-P. (2012). The area under the disease progress
stairs: calculation, advantage, and application. *Phytopathology* 102,
381-389.

## See also

[`audpc_by_unit()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc_by_unit.md)
to compute them for every plot of a trial.

## Examples

``` r
days <- c(7, 14, 21, 28, 35)
sev <- c(1, 4, 12, 30, 52)
audpc(sev, days)
#> [1] 507.5
audps(sev, days)
#> [1] 693
audpc(sev, days, type = "relative")
#> [1] 0.18125
audps(rbind(plot1 = sev, plot2 = sev / 2), days, type = "standardized")
#> plot1 plot2 
#>  19.8   9.9 
```
