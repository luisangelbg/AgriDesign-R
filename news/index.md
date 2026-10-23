# Changelog

## AgriDesign 1.0.0

First public release.

### Analysis

- [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
  analyses 16 designs (completely randomised, randomised complete blocks
  and generalised blocks, Latin, Graeco-Latin and Youden squares,
  incomplete blocks, augmented, factorials, split-, strip- and
  split-split-plots, nested) with covariates and transformations; every
  term is tested against its own error stratum.
- [`ls_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/ls_means.md),
  [`compare_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/compare_means.md)
  with eleven procedures and compact letters, simple effects with
  Satterthwaite degrees of freedom,
  [`trend()`](https://luisangelbg.github.io/AgriDesign-R/reference/trend.md),
  [`test_contrasts()`](https://luisangelbg.github.io/AgriDesign-R/reference/test_contrasts.md)
  and
  [`waller_duncan()`](https://luisangelbg.github.io/AgriDesign-R/reference/waller_duncan.md).
- [`check_assumptions()`](https://luisangelbg.github.io/AgriDesign-R/reference/check_assumptions.md),
  [`transformations()`](https://luisangelbg.github.io/AgriDesign-R/reference/transformations.md),
  [`box_cox()`](https://luisangelbg.github.io/AgriDesign-R/reference/box_cox.md),
  [`nonparametric()`](https://luisangelbg.github.io/AgriDesign-R/reference/nonparametric.md),
  [`durbin_test()`](https://luisangelbg.github.io/AgriDesign-R/reference/durbin_test.md),
  [`waerden_test()`](https://luisangelbg.github.io/AgriDesign-R/reference/waerden_test.md)
  and
  [`median_test()`](https://luisangelbg.github.io/AgriDesign-R/reference/median_test.md).
- [`recover_interblock()`](https://luisangelbg.github.io/AgriDesign-R/reference/recover_interblock.md):
  combined intra- and inter-block analysis by REML.
- [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and
  [`anova()`](https://rdrr.io/r/stats/anova.html) methods for fitted
  designs, and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for
  mean comparisons.

### Planning

- [`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)
  randomises 16 designs (uniform Latin squares by the Jacobson-Matthews
  chain; alpha and cyclic designs searched for the highest efficiency),
  [`field_book()`](https://luisangelbg.github.io/AgriDesign-R/reference/field_book.md),
  [`design_efficiency()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_efficiency.md)
  and
  [`power_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/power_design.md).
  `designs("layout")` lists the designs with the arguments each one
  needs, and every layout reports the
  [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
  call that analyses it.

### Plant breeding and pathology

- [`ammi()`](https://luisangelbg.github.io/AgriDesign-R/reference/ammi.md),
  [`stability()`](https://luisangelbg.github.io/AgriDesign-R/reference/stability.md),
  [`line_tester()`](https://luisangelbg.github.io/AgriDesign-R/reference/line_tester.md),
  [`diallel()`](https://luisangelbg.github.io/AgriDesign-R/reference/diallel.md),
  [`nc_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/nc_design.md),
  [`audpc()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc.md),
  [`audps()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc.md),
  [`audpc_by_unit()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc_by_unit.md),
  [`smith_index()`](https://luisangelbg.github.io/AgriDesign-R/reference/smith_index.md)
  and
  [`hatheway_size()`](https://luisangelbg.github.io/AgriDesign-R/reference/hatheway_size.md).

### Figures and report

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  every result,
  [`plot_interaction()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_interaction.md),
  [`plot_cell_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_cell_means.md),
  [`plot_differences()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_differences.md),
  [`plot_ss_partition()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_ss_partition.md),
  [`plot_field()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_field.md),
  [`save_figure()`](https://luisangelbg.github.io/AgriDesign-R/reference/save_figure.md)
  and
  [`design_report()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_report.md).
