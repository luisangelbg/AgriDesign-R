## Submission

This is the first submission of AgriDesign.

## Test environments

* Windows 11, R 4.6.1 (local)
* GitHub Actions: Windows (release), macOS (release), Ubuntu (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewers

* The package uses base R only (graphics, grDevices, stats, utils); knitr and
  rmarkdown are suggested for the vignettes.
* Functions that randomise layouts take a `seed` argument and restore the
  random state of the session afterwards.
* Plotting functions restore the graphical parameters they change.
* Examples run in less than 5 seconds each; tests in less than one minute.
* The Description uses British spelling (`Language: en-GB`). Words that may be
  flagged as misspelled are surnames of the authors of the methods (Dunnett,
  Eberhart, Gollob, Griffing, Piepho, Satterthwaite, Simko, Youden) or
  standard terms of experimental design (Graeco, diallel).
