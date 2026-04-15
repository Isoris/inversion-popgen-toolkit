# exactplot

**Export ggplot2 figures where the dimensions you specify are the _plot area_, not the full figure.**

## The Problem

When you call `ggsave("fig.pdf", width = 3, height = 4, units = "cm")`, that 3×4 cm is the
**entire device** — including axis labels, titles, legends, and margins. The actual plotting
area (the panel where your data lives) ends up being some unpredictable smaller size.

This means:
- You can never precisely control how big your data region is
- Two figures with different axis labels will have different panel sizes
- You can't guarantee "1 data unit = X cm" for scale-accurate plots
- Journal figure specs that refer to the plot area are impossible to meet exactly

## The Solution

```r
library(exactplot)

p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
  labs(title = "Weight vs MPG", x = "Weight (1000 lbs)", y = "Miles per gallon")

# Panel area will be EXACTLY 3 cm × 4 cm. Period.
ggsave_exact("figure1.pdf", p, panel_width = 3, panel_height = 4, units = "cm")
```

## Installation

```r
# From the tarball:
install.packages("exactplot_0.1.0.tar.gz", repos = NULL, type = "source")

# Or directly from the directory:
devtools::install("path/to/exactplot")
```

## Usage

### 1. Basic: Exact panel dimensions

```r
# 3 cm × 4 cm panel — axis text, title, legend are OUTSIDE this
ggsave_exact("fig.pdf", p, panel_width = 3, panel_height = 4, units = "cm")

# Same in inches
ggsave_exact("fig.png", p, panel_width = 1.5, panel_height = 2, units = "in", dpi = 300)

# Millimetres
ggsave_exact("fig.tiff", p, panel_width = 80, panel_height = 60, units = "mm", dpi = 600)
```

### 2. Data-coordinate scaling

"I want 1 unit on the x-axis to equal exactly 0.5 cm"

```r
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()

# Derive panel size from data range × scale factor
dims <- scale_to_data(p, x_scale = 0.5, y_scale = 0.15, units = "cm")

ggsave_exact("fig_scaled.pdf", p,
  panel_width  = dims$panel_width,
  panel_height = dims$panel_height,
  units = "cm"
)
# Now wt range (1.5–5.4 ≈ 3.9 units) × 0.5 = ~1.95 cm wide
# and mpg range (10–34 ≈ 24 units) × 0.15 = ~3.6 cm tall
```

### 3. Mixed: fix one axis, scale the other

```r
# Fixed height of 4 cm, but width follows data scale
dims <- scale_to_data(p, x_scale = 1, panel_height = 4, units = "cm")
ggsave_exact("fig_mixed.pdf", p,
  panel_width = dims$panel_width, panel_height = dims$panel_height, units = "cm")
```

### 4. Inspect dimensions without saving

```r
dims <- exact_dims(p, panel_width = 3, panel_height = 4, units = "cm")
dims$device_width   # total device width in inches (what ggsave would need)
dims$device_height
dims$chrome$left    # space eaten by y-axis label + tick labels
dims$chrome$right   # space on the right (legend if present)
dims$chrome$top     # title
dims$chrome$bottom  # x-axis label + tick labels
```

### 5. Verbose mode

```r
ggsave_exact("fig.pdf", p, panel_width = 3, panel_height = 4,
  units = "cm", verbose = TRUE)
# exactplot: panel target = 1.181 x 1.575 in | device = 2.019 x 2.512 in | iterations = 2
#   chrome: L=0.486 R=0.352 T=0.446 B=0.491 in
```

### 6. Unit conversion helper

```r
unit_convert(2.54, "cm", "in")  # 1
unit_convert(1, "in", "pt")     # 72
unit_convert(72, "pt", "mm")    # 25.4
```

## How It Works

1. **Trial render**: Opens a null graphics device and builds the ggplot into a gtable
2. **Measure chrome**: Sums the widths/heights of all grob rows/columns outside the panel cells
3. **Inflate**: device_size = panel_target + chrome
4. **Iterate**: Re-renders at the new device size (because chrome can shift slightly when the device changes) until convergence (default tolerance: 0.005 inches ≈ 0.13 mm)
5. **Save**: Delegates to `ggplot2::ggsave()` with the computed device dimensions

## Supports

- All ggplot2 plot types (scatter, bar, line, boxplot, etc.)
- Faceted plots (measures the bounding box of all panels)
- All output formats supported by ggsave (PDF, PNG, TIFF, SVG, etc.)
- Custom themes with varying margin sizes
- Plots with or without legends, titles, subtitles, captions

## Limitations

- **Facets with `scales = "free"`**: The panel_width/height applies to the bounding box of all panels combined, not individual facet panels.
- **Interactive/dynamic plots**: This is for static export only.
- **coord_polar / coord_sf**: May not work perfectly with non-Cartesian coordinate systems.

## License

MIT
