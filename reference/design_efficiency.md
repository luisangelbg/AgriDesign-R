# Efficiency factor of a block design

The efficiency factor of a block design is the harmonic mean of its
canonical efficiency factors: the average variance of a treatment
difference in a complete block design with the same replication, divided
by the average variance in this design (intra-block analysis, same error
variance). It is 1 for complete blocks and the target of incomplete
block designs. For resolvable designs the upper bound of Patterson and
Williams (1976) is also reported.

## Usage

``` r
design_efficiency(x, treatment = "Treatment", block = NULL)
```

## Arguments

- x:

  An `"ad_layout"` object from
  [`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md),
  or a data frame.

- treatment, block:

  Column names of the treatment and of the (incomplete) block when `x`
  is a data frame. For a layout the blocks are the incomplete blocks
  when there are any.

## Value

A list with the efficiency factor `E`, the canonical efficiency factors,
whether the design is connected, and `upper_bound` for equireplicated
resolvable designs.

## Examples

``` r
lay <- make_design("alpha", trt = 12, r = 2, k = 3, seed = 1)
design_efficiency(lay)$E
#> [1] 0.5945946
```
