# Save a figure to a file

Opens the graphics device that matches the extension of `file` (`.png`,
`.tiff`, `.jpg`, `.pdf`, `.svg` or `.eps`), evaluates the plotting code
and closes the device. Raster formats carry the resolution in their
header, as journals ask.

## Usage

``` r
save_figure(
  file,
  expr,
  width = 16,
  height = 11,
  units = c("cm", "mm", "in"),
  res = 600,
  pointsize = 10
)
```

## Arguments

- file:

  Output file; its extension chooses the format.

- expr:

  Plotting code, for instance `plot(cmp)`; it is evaluated with the
  device open.

- width, height:

  Size of the figure.

- units:

  `"cm"`, `"mm"` or `"in"`.

- res:

  Resolution in dots per inch of raster formats.

- pointsize:

  Base size of the text, in points.

## Value

The path of the file, invisibly.

## Examples

``` r
fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
f <- tempfile(fileext = ".png")
save_figure(f, plot(compare_means(fit, "Nitrogen_kg_ha")), width = 12, height = 9)
file.exists(f)
#> [1] TRUE
```
