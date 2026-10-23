# AgriDesign 1.0.0

First public release.

## Analysis

* `design_anova()` analyses 16 designs (completely randomised, randomised
  complete blocks and generalised blocks, Latin, Graeco-Latin and Youden
  squares, incomplete blocks, augmented, factorials, split-, strip- and
  split-split-plots, nested) with covariates and transformations; every term
  is tested against its own error stratum.
* `ls_means()`, `compare_means()` with eleven procedures and compact letters,
  simple effects with Satterthwaite degrees of freedom, `trend()`,
  `test_contrasts()` and `waller_duncan()`.
* `check_assumptions()`, `transformations()`, `box_cox()`, `nonparametric()`,
  `durbin_test()`, `waerden_test()` and `median_test()`.
* `recover_interblock()`: combined intra- and inter-block analysis by REML.
* `residuals()`, `fitted()`, `nobs()` and `anova()` methods for fitted designs, and
  `as.data.frame()` for mean comparisons.

## Planning

* `make_design()` randomises 16 designs (uniform Latin squares by the
  Jacobson-Matthews chain; alpha and cyclic designs searched for the highest
  efficiency), `field_book()`, `design_efficiency()` and `power_design()`.
  `designs("layout")` lists the designs with the arguments each one needs,
  and every layout reports the `design_anova()` call that analyses it.

## Plant breeding and pathology

* `ammi()`, `stability()`, `line_tester()`, `diallel()`, `nc_design()`,
  `audpc()`, `audps()`, `audpc_by_unit()`, `smith_index()` and
  `hatheway_size()`.

## Figures and report

* `plot()` methods for every result, `plot_interaction()`,
  `plot_cell_means()`, `plot_differences()`, `plot_ss_partition()`,
  `plot_field()`, `save_figure()` and `design_report()`.
