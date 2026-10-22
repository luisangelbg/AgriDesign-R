#' Maize yield under four nitrogen rates (RCBD)
#'
#' Randomised complete block design with four nitrogen rates in four blocks.
#'
#' @format A data frame with 16 rows and 3 columns:
#' \describe{
#'   \item{Block}{Block (1-4).}
#'   \item{Nitrogen_kg_ha}{Nitrogen rate, kg/ha (0, 60, 120, 180).}
#'   \item{Yield_t_ha}{Grain yield, t/ha.}
#' }
#' @examples
#' design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
"maize_nitrogen"

#' Yield of five bean varieties (CRD)
#'
#' Completely randomised design with five common bean varieties.
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{Plot}{Plot number.}
#'   \item{Variety}{Bean variety.}
#'   \item{Yield_kg_plot}{Yield, kg per plot.}
#' }
#' @examples
#' design_anova(bean_varieties, "Yield_kg_plot", "Variety")
"bean_varieties"

#' Disease severity scores under four fungicides (RCBD)
#'
#' Ordinal severity scores (0-9) of a foliar disease under four fungicide
#' treatments in five blocks; suited to the rank-based (Friedman) analysis.
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{Block}{Block (1-5).}
#'   \item{Fungicide}{Fungicide treatment, including an untreated control.}
#'   \item{Severity_score}{Disease severity score.}
#' }
"disease_scores"

#' Tomato varieties by fertiliser (factorial in RCBD)
#'
#' Three tomato varieties by three fertiliser levels in three blocks.
#'
#' @format A data frame with 27 rows and 5 columns:
#' \describe{
#'   \item{Block}{Block (1-3).}
#'   \item{Variety}{Tomato variety.}
#'   \item{Fertilizer}{Fertiliser level (F0, F1, F2).}
#'   \item{Yield_kg_plant}{Fruit yield, kg per plant.}
#'   \item{Height_cm}{Plant height, cm.}
#' }
#' @examples
#' design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"),
#'              block = "Block")
"tomato_factorial"

#' Wheat yield in a 4 x 4 Latin square
#'
#' @format A data frame with 16 rows and 4 columns:
#' \describe{
#'   \item{Row}{Row of the field (1-4).}
#'   \item{Column}{Column of the field (1-4).}
#'   \item{Treatment}{Treatment (A-D).}
#'   \item{Yield_t_ha}{Grain yield, t/ha.}
#' }
#' @examples
#' design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column")
"wheat_latin"

#' Irrigation systems by varieties (split-plot in RCBD)
#'
#' Two irrigation systems on main plots and four varieties on sub-plots, in
#' three blocks.
#'
#' @format A data frame with 24 rows and 4 columns:
#' \describe{
#'   \item{Block}{Block (1-3).}
#'   \item{Irrigation}{Irrigation system (main-plot factor).}
#'   \item{Variety}{Variety (sub-plot factor).}
#'   \item{Yield_t_ha}{Yield, t/ha.}
#' }
#' @examples
#' design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"),
#'              block = "Block", design = "split_rcbd")
"split_irrigation"

#' Maize nitrogen trial in wide format
#'
#' The same trial as [maize_nitrogen], with one column per nitrogen rate.
#' It shows the wide layout that has to be reshaped to long format (one row
#' per plot) before the analysis.
#'
#' @format A data frame with 4 rows and 5 columns: `Block` and the yields
#'   (t/ha) at `N0`, `N60`, `N120` and `N180`.
"maize_wide"
