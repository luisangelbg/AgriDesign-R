# Package index

## Analysis of variance

Fit the model of a design and get its least-squares means.

- [`design_anova()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_anova.md)
  : Analysis of variance of an experimental design
- [`designs()`](https://luisangelbg.github.io/AgriDesign-R/reference/designs.md)
  : List the designs that can be analysed or randomised
- [`ls_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/ls_means.md)
  : Least-squares means of a design
- [`recover_interblock()`](https://luisangelbg.github.io/AgriDesign-R/reference/recover_interblock.md)
  : Recover inter-block information (REML)
- [`residuals(`*`<ad_anova>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/ad_methods.md)
  [`fitted(`*`<ad_anova>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/ad_methods.md)
  [`nobs(`*`<ad_anova>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/ad_methods.md)
  [`anova(`*`<ad_anova>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/ad_methods.md)
  [`as.data.frame(`*`<ad_comparison>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/ad_methods.md)
  : Standard methods for fitted designs

## Assumptions

- [`check_assumptions()`](https://luisangelbg.github.io/AgriDesign-R/reference/check_assumptions.md)
  : Check the assumptions of the analysis of variance
- [`transformations()`](https://luisangelbg.github.io/AgriDesign-R/reference/transformations.md)
  : Compare transformations of the response
- [`box_cox()`](https://luisangelbg.github.io/AgriDesign-R/reference/box_cox.md)
  : Box-Cox profile likelihood of a design

## Comparing treatments

- [`compare_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/compare_means.md)
  : Compare treatment means
- [`waller_duncan()`](https://luisangelbg.github.io/AgriDesign-R/reference/waller_duncan.md)
  : Waller-Duncan k-ratio t test
- [`waller_critical()`](https://luisangelbg.github.io/AgriDesign-R/reference/waller_critical.md)
  : Critical t of the Waller-Duncan k-ratio test
- [`trend()`](https://luisangelbg.github.io/AgriDesign-R/reference/trend.md)
  : Trend analysis of a quantitative factor
- [`test_contrasts()`](https://luisangelbg.github.io/AgriDesign-R/reference/test_contrasts.md)
  : Test contrasts among the levels of a factor

## Rank-based tests

- [`nonparametric()`](https://luisangelbg.github.io/AgriDesign-R/reference/nonparametric.md)
  : Rank-based and unequal-variance alternatives to the ANOVA
- [`durbin_test()`](https://luisangelbg.github.io/AgriDesign-R/reference/durbin_test.md)
  : Durbin test for balanced incomplete blocks
- [`waerden_test()`](https://luisangelbg.github.io/AgriDesign-R/reference/waerden_test.md)
  : van der Waerden normal-scores test
- [`median_test()`](https://luisangelbg.github.io/AgriDesign-R/reference/median_test.md)
  : Mood's median test

## Planning experiments

- [`make_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/make_design.md)
  : Randomise a field experiment
- [`field_book()`](https://luisangelbg.github.io/AgriDesign-R/reference/field_book.md)
  : Field book of a randomised layout
- [`design_efficiency()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_efficiency.md)
  : Efficiency factor of a block design
- [`power_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/power_design.md)
  : Replicates and power for a difference between two treatments

## Plant breeding

- [`ammi()`](https://luisangelbg.github.io/AgriDesign-R/reference/ammi.md)
  : AMMI analysis of a multi-environment trial
- [`stability()`](https://luisangelbg.github.io/AgriDesign-R/reference/stability.md)
  : Stability statistics of a multi-environment trial
- [`line_tester()`](https://luisangelbg.github.io/AgriDesign-R/reference/line_tester.md)
  : Line x tester analysis
- [`diallel()`](https://luisangelbg.github.io/AgriDesign-R/reference/diallel.md)
  : Diallel analysis by Griffing's methods
- [`nc_design()`](https://luisangelbg.github.io/AgriDesign-R/reference/nc_design.md)
  : North Carolina designs I, II and III

## Plant pathology and plot size

- [`audpc()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc.md)
  [`audps()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc.md)
  : Area under the disease progress curve and stairs
- [`audpc_by_unit()`](https://luisangelbg.github.io/AgriDesign-R/reference/audpc_by_unit.md)
  : AUDPC and AUDPS of every plot of a trial
- [`smith_index()`](https://luisangelbg.github.io/AgriDesign-R/reference/smith_index.md)
  : Smith's index of soil heterogeneity
- [`hatheway_size()`](https://luisangelbg.github.io/AgriDesign-R/reference/hatheway_size.md)
  : Plot size needed to detect a difference (Hatheway)

## Figures and report

- [`plot(`*`<ad_ammi>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/plot.ad_ammi.md)
  : Biplots of an AMMI analysis
- [`plot(`*`<ad_anova>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/plot.ad_anova.md)
  : Residual plots of a design
- [`plot(`*`<ad_comparison>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/plot.ad_comparison.md)
  : Plot treatment means with their letters
- [`plot(`*`<ad_layout>`*`)`](https://luisangelbg.github.io/AgriDesign-R/reference/plot.ad_layout.md)
  : Plot a randomised field layout
- [`plot_cell_means()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_cell_means.md)
  : Heat map of cell means
- [`plot_differences()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_differences.md)
  : Plot the pairwise differences with their intervals
- [`plot_field()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_field.md)
  : Field map of the response or of the residuals
- [`plot_interaction()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_interaction.md)
  : Interaction plot of two treatment factors
- [`plot_ss_partition()`](https://luisangelbg.github.io/AgriDesign-R/reference/plot_ss_partition.md)
  : Partition of the total sum of squares
- [`save_figure()`](https://luisangelbg.github.io/AgriDesign-R/reference/save_figure.md)
  : Save a figure to a file
- [`agri_colors()`](https://luisangelbg.github.io/AgriDesign-R/reference/agri_colors.md)
  : Colours of the AgriDesign figures
- [`design_report()`](https://luisangelbg.github.io/AgriDesign-R/reference/design_report.md)
  : Write an HTML report of an analysis

## Example data

- [`maize_nitrogen`](https://luisangelbg.github.io/AgriDesign-R/reference/maize_nitrogen.md)
  : Maize yield under four nitrogen rates (RCBD)
- [`maize_wide`](https://luisangelbg.github.io/AgriDesign-R/reference/maize_wide.md)
  : Maize nitrogen trial in wide format
- [`bean_varieties`](https://luisangelbg.github.io/AgriDesign-R/reference/bean_varieties.md)
  : Yield of five bean varieties (CRD)
- [`disease_scores`](https://luisangelbg.github.io/AgriDesign-R/reference/disease_scores.md)
  : Disease severity scores under four fungicides (RCBD)
- [`tomato_factorial`](https://luisangelbg.github.io/AgriDesign-R/reference/tomato_factorial.md)
  : Tomato varieties by fertiliser (factorial in RCBD)
- [`wheat_latin`](https://luisangelbg.github.io/AgriDesign-R/reference/wheat_latin.md)
  : Wheat yield in a 4 x 4 Latin square
- [`split_irrigation`](https://luisangelbg.github.io/AgriDesign-R/reference/split_irrigation.md)
  : Irrigation systems by varieties (split-plot in RCBD)
- [`bib_sorghum`](https://luisangelbg.github.io/AgriDesign-R/reference/bib_sorghum.md)
  : Sorghum hybrids in a balanced incomplete block design (simulated)
- [`blight_progress`](https://luisangelbg.github.io/AgriDesign-R/reference/blight_progress.md)
  : Late blight progress under four fungicide programmes (simulated)
- [`uniformity_trial`](https://luisangelbg.github.io/AgriDesign-R/reference/uniformity_trial.md)
  : Uniformity trial of 24 x 16 basic units (simulated)
- [`maize_met`](https://luisangelbg.github.io/AgriDesign-R/reference/maize_met.md)
  : Maize hybrids in six environments (multi-environment trial)
- [`maize_lxt`](https://luisangelbg.github.io/AgriDesign-R/reference/maize_lxt.md)
  : Line x tester trial in maize
- [`maize_diallel`](https://luisangelbg.github.io/AgriDesign-R/reference/maize_diallel.md)
  : Full diallel among six maize lines
- [`nc_designs`](https://luisangelbg.github.io/AgriDesign-R/reference/nc_designs.md)
  : North Carolina mating designs

## Package

- [`AgriDesign`](https://luisangelbg.github.io/AgriDesign-R/reference/AgriDesign-package.md)
  [`AgriDesign-package`](https://luisangelbg.github.io/AgriDesign-R/reference/AgriDesign-package.md)
  : AgriDesign: design and analysis of agricultural experiments
