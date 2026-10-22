# Stage 6 validation: HTML reports of the example data sets. Checks that the
# numbers printed are those of the analysis, that the page is well formed and
# self-contained, and that the inline SVG figures do not share element ids.
suppressMessages(devtools::load_all(quiet = TRUE))
out <- "validation/reports"; dir.create(out, showWarnings = FALSE)
ok <- 0; bad <- 0
chk <- function(what, pass, info = "") {
  pass <- isTRUE(all(pass)); if (pass) ok <<- ok + 1 else bad <<- bad + 1
  cat(sprintf("%-4s %s %s\n", if (pass) "OK" else "FAIL", what, info))
}
has <- function(html, s) grepl(s, html, fixed = TRUE)

fits <- list(
  maize = list(design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block"), units = "t/ha"),
  beans = list(design_anova(bean_varieties, "Yield_kg_plot", "Variety"), units = "kg per plot"),
  latin = list(design_anova(wheat_latin, "Yield_t_ha", "Treatment", row = "Row", col = "Column"), units = "t/ha"),
  tomato = list(design_anova(tomato_factorial, "Yield_kg_plant", c("Variety", "Fertilizer"), block = "Block", transform = "log"), units = "kg/plant"),
  split = list(design_anova(split_irrigation, "Yield_t_ha", c("Irrigation", "Variety"), block = "Block", design = "split_rcbd"), units = "t/ha"))
html <- list()
for (n in names(fits)) {
  f <- fits[[n]][[1]]
  t0 <- Sys.time()
  html[[n]] <- design_report(f, file.path(out, paste0(n, ".html")), units = fits[[n]]$units, data = TRUE,
                             authors = "L. A. Barrera-Guzman", trial = paste("Example data:", n))
  cat(sprintf("INFO %-7s %3d KB, %d tables, %d figures, %.1f s\n", n, round(nchar(html[[n]], "bytes") / 1024),
              attr(html[[n]], "tables"), attr(html[[n]], "figures"), as.numeric(Sys.time() - t0, units = "secs")))
}
h <- html$maize
chk("maize: F of the ANOVA printed (437.267)", has(h, ">437.267<"))
chk("maize: CV 2.1 % in the caption", has(h, "CV = 2.1 %"))
chk("maize: LS means 7.183 a ... 4.215 d", has(h, ">7.183<") && has(h, ">4.215<") && has(h, ">a<") && has(h, ">d<"))
chk("maize: HSD 0.283 in the caption", has(h, "Minimum significant difference = 0.283"))
chk("maize: trend SS 19.7707 and 1.7756", has(h, ">19.7707<") && has(h, ">1.7756<"))
chk("maize: optimum 179.5", has(h, "maximum at Nitrogen_kg_ha = 179.5"))
chk("maize: grade A box", has(h, "<span class=\"grade\">A</span>"))
chk("maize: methods name the design, model and test", has(h, "randomised complete block design") && has(h, "Yield_t_ha ~ Block + Nitrogen_kg_ha") &&
      has(h, "Tukey HSD"))
chk("maize: narrative F<sub>3,9</sub>", has(h, "F<sub>3,9</sub> = 437.27"))
h <- html$split
chk("split: Irrigation tested against Error a", has(h, "Error a (main plot)") && has(h, ">10460.641<"))
chk("split: simple effects table with pooled error", has(h, "pooled error [(b-1)Eb + Ea]/b"))
chk("split: methods mention the error strata", has(h, "testing each term against its own error stratum"))
h <- html$tomato
chk("tomato: transformation stated and back-transformed column", has(h, "Natural logarithm ln(y)") && has(h, ">Back-transformed<"))
h <- html$latin
chk("latin: rows x columns in the methods", has(h, "4 rows &times; 4 columns"))
# structure of every page
for (n in names(html)) {
  h <- html[[n]]
  opened <- lengths(regmatches(h, gregexpr("<table>", h, fixed = TRUE))); closed <- lengths(regmatches(h, gregexpr("</table>", h, fixed = TRUE)))
  ids <- regmatches(h, gregexpr("id=\"[^\"]+\"", h))[[1]]
  refs <- unique(sub("^(href=\"#|url\\(#)", "", regmatches(h, gregexpr("(href=\"#|url\\(#)[^\")]+", h))[[1]]))
  ids_clean <- sub("^id=\"", "", sub("\"$", "", ids))
  chk(sprintf("%s: tables closed, ids unique, every reference resolves", n),
      opened == closed && !anyDuplicated(ids) && all(refs %in% ids_clean),
      sprintf("(%d tables, %d ids)", opened, length(ids)))
  chk(sprintf("%s: self-contained (no external scripts, styles or images)", n),
      !grepl("<script|<link|src=\"http", h))
  chk(sprintf("%s: no entity escaped twice", n), !grepl("&amp;[a-z]+;", h))
}
chk("Source of the report is ASCII", !any(grepl("[^\\x01-\\x7F]", readLines("R/report.R"), perl = TRUE)))
cat(sprintf("\n%d checks passed, %d failed\n", ok, bad))
