# Builds data/*.rda from the example CSV files of the AgriDesign web application.
# Run from the package root:  Rscript data-raw/make-data.R
src <- "data-raw/csv"
read <- function(f) utils::read.csv(file.path(src, f), stringsAsFactors = FALSE, check.names = FALSE)

maize_nitrogen   <- read("rcbd_maize_nitrogen.csv")
bean_varieties   <- read("crd_bean_varieties.csv")
disease_scores   <- read("disease_scores_rcbd.csv")
tomato_factorial <- read("factorial_tomato_rcbd.csv")
wheat_latin      <- read("latin_square_wheat.csv")
split_irrigation <- read("split_plot_irrigation.csv")
maize_wide       <- read("wide_format_example.csv")

for (nm in c("maize_nitrogen", "bean_varieties", "disease_scores", "tomato_factorial",
             "wheat_latin", "split_irrigation", "maize_wide"))
  save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz", version = 2)
