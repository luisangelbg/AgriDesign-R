# Validation

Every statistic of AgriDesign is checked against an independent
implementation: functions of base R, established R packages, worked
examples of textbooks and software manuals, and closed-form results. The
checks are scripts in the `validation/` folder of the source repository;
each prints one line per comparison and the number of checks passed. At
the time of writing all **540** checks pass. The package itself depends
on no other package: the reference packages are used only by these
scripts.

| Script | Checks | Compared with |
|----|---:|----|
| `v01-anova.R` | 29 | [`aov()`](https://rdrr.io/r/stats/aov.html) with `Error()` strata, car (Type III), emmeans, agricolae (`sp.plot`) |
| `v02-assumptions-nonpar.R` | 38 | `shapiro.test`, nortest, car (Levene, Durbin-Watson), MASS (`boxcox`), `friedman.test`, `kruskal.test`, `oneway.test`, PMCMRplus, ARTool, rcompanion |
| `v03-means.R` | 47 | `TukeyHSD`, `pairwise.t.test`, agricolae, mvtnorm, PMCMRplus, emmeans, lmerTest |
| `v04-figures.R` | 17 | the numbers of each analysis against what is drawn; file formats and resolution |
| `v05-generator.R` | 31 | structure over many seeds, uniformity tests, `power.t.test` |
| `v06-report.R` | 30 | the figures of the analysis in the HTML, well-formed and self-contained pages |
| `v07a-designs.R` | 63 | closed forms of efficiency factors, agricolae designs, lme4 + emmeans (REML) |
| `v07b-breeding.R` | 199 | agricolae, `lm`, `svd`, lme4, exact enumeration of F2 populations and of rank distributions |
| `v07c-tests-pathology.R` | 86 | a SAS manual example, an independent double integral, PMCMRplus, `chisq.test`, `fisher.test`, agricolae |

## Differences with other software

When two programs disagree, the validation looks for a third,
independent answer. These differences were found, and in each case the
reason is stated in the validation script:

- **Dunnett’s critical value.** With its default tolerance, multcomp
  computes the multivariate t by quasi-Monte Carlo with an error of
  about 1e-3 (2.7251 instead of 2.7273 for four comparisons and 15 error
  degrees of freedom). AgriDesign integrates it with Gauss-Legendre
  quadrature, using the exact correlation of every pair, and agrees with
  mvtnorm at high precision.
- **Efficiency of alpha designs.** `agricolae::design.alpha()` prints
  the upper bound of Patterson and Williams as the efficiency factor of
  the design it builds. For 12 treatments in blocks of 3 with 2
  replicates it prints 0.6471, while its layout has an efficiency of
  0.5946.
- **Latin squares.** Permuting the rows, columns and symbols of a cyclic
  square reaches only 432 of the 576 Latin squares of order 4.
  AgriDesign draws them uniformly with the Jacobson-Matthews chain.
- **North Carolina design III.** The dominance variance equals the males
  x testers component; doubling it (as agricolae does) overestimates it
  by a factor of two, as the exact enumeration of an F2 population
  shows.
- **Line x tester.** The average covariance of half-sibs needs the
  factor (l + t - 2) of the weighted mean square, the covariance of
  full-sibs keeps every term of the formula of Singh and Chaudhary, and
  the dominance variance is 4 sigma2_SCA / (1 + F)^2.
- **Huehn’s S(1).** Its variance includes the factor (K^2 - 1) / K^2, as
  the exact null distribution confirms.
- **Waller-Duncan.** The critical t agrees with the SAS manual example
  and with an independent double integral; agricolae’s approximation
  drifts when the treatments have one degree of freedom or when F is
  large.
- **Median test.** `agricolae::Median.test(correct = TRUE)` does not
  apply the continuity correction, because the argument reaches
  [`chisq.test()`](https://rdrr.io/r/stats/chisq.test.html) in the wrong
  position.
- **van der Waerden test.** Rounding the ranks to three decimals before
  the normal scores changes the statistic in the fourth significant
  digit.
- **Power.** A shifted central t underestimates the power with few error
  degrees of freedom (0.21 instead of 0.245 for a CV of 15 %, a
  difference of 15 % and four replicates); AgriDesign uses the
  noncentral t.
