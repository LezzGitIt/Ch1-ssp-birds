## PhD birds in silvopastoral landscapes ##
# Land-use context for the habitat-by-time interaction: forest area within 5 km of the Meta point counts from 2000 to 2024, from Hansen Global Forest Change (annual loss) and the Colombia woody canopy cover and height maps (survey years), to show whether the landscape around the farms lost or gained forest over the study period.

### The 5 km buffers around every Meta point are dissolved into one polygon before summing, so pixels where neighbouring buffers overlap are counted once.
### Forest follows Hansen's definition in both sources: canopy cover of at least 30% in vegetation at least 5 m tall.
## Hansen: forest in 2000 (treecover2000 at or above the threshold) plus the 2000–2012 gain, minus the cumulative stand-replacing loss up to each year; the gain layer is a single undated 2000–2012 flag, so the series is drawn as a straight line from 2000 to 2012 and only falls with later loss from 2013 on.
## Woody canopy maps (Potapov 2025): pixels with cover and woody canopy height both at or above the thresholds in each survey raster year, so this series can show gains as well as losses.

# Load libraries ----
library(tidyverse)
library(cowplot)
ggplot2::theme_set(cowplot::theme_cowplot())

# Configuration ----
Buffer_m <- 5000
Forest_cover_pct <- 30
Forest_height_m <- 5

## Hansen Global Forest Change v1.12 (2001–2024), 30 m; tile 10N_080W (0–10°N, 80–70°W) holds every Meta point
Hansen_version <- "GFC-2024-v1.12"
Hansen_tile <- "10N_080W"
Hansen_layers <- c("treecover2000", "lossyear", "gain")
Hansen_dir <- "../Geospatial_data/Environmental/Hansen_GFC_v1.12"
Clip_buffer_m <- 20000

## Colombia woody canopy cover (WVCC, %) and height (WVCH, decimetres) maps for the survey years
National_dir <- "../Geospatial_data/Environmental"
Raster_years <- c(2013, 2014, 2016, 2017, 2019, 2022, 2024)

## Colombia Woody Vegetation Structure Change (WVSC, Zenodo record 18154841): one map per year, 2010–2024, classing each pixel's change in woody vegetation (trees and shrubs ≥ 2 m) from the previous year; kept as clips to a dissolved 50 km buffer around the Meta points. The product documentation warns that 2010–2014 are less precise and that gradual increases may be inconsistent before vs after 2015.
Wvsc_years <- 2010:2024
Wvsc_dir <- file.path(National_dir, "WVSC")
Wvsc_clip_buffer_m <- 50000
Wvsc_classes <- c("1" = "Removal", "2" = "Degradation", "3" = "Establishment", "4" = "Enhancement")

## Point locations and survey dates from the wrangling repo; buffers are drawn in UTM zone 18N
Wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"
Utm_crs <- "EPSG:32618"

## Series labels in the table and legend
Hansen_label <- "Hansen et al. (2013)"
Potapov_label <- "Potapov (2025)"

## Outputs
Forest_csv <- "Derived/Excels/Forest_area_meta_5km.csv"
Woody_change_csv <- "Derived/Excels/Woody_change_meta_5km.csv"
Fig_dir <- "Figures/Land_use_context"

# Functions ----
## Crop a raster to a polygon and set cells outside it to NA
read_in_area <- function(raster_path, area) {
  terra::rast(raster_path) %>%
    terra::crop(area, snap = "out") %>%
    terra::mask(area)
}

## Bar panel of annual change by class, from a table with Year, Class and Signed_ha (losses negative); y axis labels show absolute areas
change_panel <- function(change, title, fill_values, x_breaks, y_limits = NULL, y_breaks = waiver()) {
  ggplot(change, aes(Year, Signed_ha, fill = Class)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = fill_values) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(labels = abs, limits = y_limits, breaks = y_breaks) +
    labs(title = title, x = "Year", y = "Hectares", fill = NULL) +
    theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(size = 12, face = "plain"))
}

## Small line plot of cumulative net change in area (Year, Net_ha) for an inset, labelled at its last point; y_limits lets several insets share one scale
net_change_inset <- function(net, title, label, x_breaks, y_limits) {
  ggplot(net, aes(Year, Net_ha)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    annotate("text", x = max(net$Year), y = last(net$Net_ha), label = label, hjust = 1, vjust = -0.6, size = 2.6) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = 0.25)) +
    labs(title = title, x = NULL, y = NULL) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(plot.background = element_rect(fill = "white", colour = "grey70"), plot.margin = margin(3, 5, 3, 3), plot.title = element_text(size = 8, face = "plain"))
}

# Buffer the Meta points ----
meta_pts <- read_csv(paste0(Wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  filter(Department == "Meta")
dissolved_buffer <- function(points, radius_m) {
  terra::vect(points, geom = c("Long", "Lat"), crs = "EPSG:4326") %>%
    terra::project(Utm_crs) %>%
    terra::buffer(radius_m) %>%
    terra::aggregate() %>%
    terra::project("EPSG:4326")
}
buffer_union <- dissolved_buffer(meta_pts, Buffer_m)

# Clip Hansen tiles ----
## Each layer is kept only within a dissolved Clip_buffer_m buffer around the Meta points (a few MB instead of ~630 MB for the full tiles); a missing clip is made from the full tile, read from Hansen_dir if present and otherwise downloaded to a temporary file
hansen_paths <- set_names(file.path(Hansen_dir, sprintf("Hansen_%s_%s_Meta_%dkm.tif", Hansen_version, Hansen_layers, Clip_buffer_m / 1000)), Hansen_layers)
dir.create(Hansen_dir, showWarnings = FALSE)
options(timeout = 3600)
clip_area <- dissolved_buffer(meta_pts, Clip_buffer_m)
walk(Hansen_layers[!file.exists(hansen_paths)], \(layer) {
  tile_name <- sprintf("Hansen_%s_%s_%s.tif", Hansen_version, layer, Hansen_tile)
  tile_path <- file.path(Hansen_dir, tile_name)
  if (!file.exists(tile_path)) {
    tile_path <- tempfile(fileext = ".tif")
    download.file(paste0("https://storage.googleapis.com/earthenginepartners-hansen/", Hansen_version, "/", tile_name), tile_path, mode = "wb")
  }
  terra::writeRaster(read_in_area(tile_path, clip_area), hansen_paths[[layer]], datatype = "INT1U", gdal = "COMPRESS=LZW")
})

# Clip WVSC maps ----
## Each missing year is downloaded to a temporary file, clipped, and the national file deleted
wvsc_paths <- set_names(file.path(Wvsc_dir, sprintf("Colombia_WVSC_%d_Meta_%dkm.tif", Wvsc_years, Wvsc_clip_buffer_m / 1000)), Wvsc_years)
dir.create(Wvsc_dir, showWarnings = FALSE)
wvsc_clip_area <- dissolved_buffer(meta_pts, Wvsc_clip_buffer_m)
walk(Wvsc_years[!file.exists(wvsc_paths)], \(year) {
  national_path <- tempfile(fileext = ".tif")
  download.file(sprintf("https://zenodo.org/records/18154841/files/Colombia_WVSC_%d.tif?download=1", year), national_path, mode = "wb")
  terra::writeRaster(read_in_area(national_path, wvsc_clip_area), wvsc_paths[[as.character(year)]], datatype = "INT1U", gdal = "COMPRESS=LZW")
  unlink(national_path)
})

# Hansen forest area by year ----
hansen <- map(hansen_paths, read_in_area, area = buffer_union)
hansen_cells <- tibble(
  Area_ha = terra::values(terra::cellSize(hansen$treecover2000, unit = "ha", mask = TRUE))[, 1],
  Treecover2000 = terra::values(hansen$treecover2000)[, 1],
  Lossyear = terra::values(hansen$lossyear)[, 1],
  Gain = terra::values(hansen$gain)[, 1]
) %>%
  filter(!is.na(Area_ha))
hansen_buffer_ha <- sum(hansen_cells$Area_ha)
forest_2000 <- hansen_cells %>% filter(Treecover2000 >= Forest_cover_pct)
hansen_gain_ha <- sum(hansen_cells$Area_ha[hansen_cells$Gain == 1])
## Gain cannot be dated within 2000–2012, so the series has no values between 2000 and 2012 (a straight line in the figure, gain spread evenly); forest in 2012 is every pixel still forest from 2000 or flagged as gain, and pixels flagged with both gain and loss by 2012 count as forest in 2012
forest_2012 <- hansen_cells %>%
  filter(Gain == 1 | (Treecover2000 >= Forest_cover_pct & (Lossyear == 0 | Lossyear > 12)))
hansen_forest <- tibble(Year = c(2000, 2012:2024)) %>%
  mutate(
    Source = Hansen_label,
    Forest_ha = map_dbl(Year, \(year) if (year == 2000) sum(forest_2000$Area_ha) else sum(forest_2012$Area_ha[forest_2012$Lossyear <= 12 | forest_2012$Lossyear > year - 2000])),
    Buffer_ha = hansen_buffer_ha
  )

# Woody canopy forest area by survey year ----
woody_forest <- map(Raster_years, \(year) {
  cover <- read_in_area(file.path(National_dir, sprintf("WVCC/Colombia_WVCC_%d_50km.tif", year)), buffer_union)
  height <- read_in_area(file.path(National_dir, sprintf("WVCH/Colombia_WVCH_%d_50km.tif", year)), buffer_union)
  cell_ha <- terra::cellSize(cover, unit = "ha", mask = TRUE)
  forest <- cover >= Forest_cover_pct & height >= Forest_height_m * 10
  tibble(
    Year = year,
    Source = Potapov_label,
    Forest_ha = terra::global(cell_ha * forest, "sum", na.rm = TRUE)[[1]],
    Buffer_ha = terra::global(cell_ha, "sum", na.rm = TRUE)[[1]]
  )
}) %>%
  list_rbind()

# Save table ----
forest_area <- bind_rows(hansen_forest, woody_forest) %>%
  mutate(Forest_pct = 100 * Forest_ha / Buffer_ha) %>%
  select(Source, Year, Forest_ha, Buffer_ha, Forest_pct) %>%
  arrange(Source, Year) %>%
  mutate(across(c(Forest_ha, Buffer_ha), \(x) round(x, 1)), Forest_pct = round(Forest_pct, 2))
write_csv(forest_area, Forest_csv)

# Plot ----
## Shaded bands mark the years each survey year-group spans in Meta
year_groups <- read_csv(paste0(Wrangling_excels, "Event_covs.csv"), show_col_types = FALSE) %>%
  semi_join(meta_pts, by = "Id_survey_no_dc") %>%
  group_by(Year_grp) %>%
  summarise(Start = min(Year) - 0.5, End = max(Year) + 0.5)
forest_plot <- ggplot(forest_area, aes(Year, Forest_pct, colour = Source)) +
  geom_rect(data = year_groups, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_line() +
  geom_point(data = filter(forest_area, Source == Potapov_label)) +
  scale_colour_manual(values = set_names(c("#1b7837", "#762a83"), c(Hansen_label, Potapov_label))) +
  labs(x = "Year", y = "Forest within 5 km (%)", colour = NULL) +
  theme(legend.position = "bottom")
dir.create(Fig_dir, showWarnings = FALSE)
ggsave(file.path(Fig_dir, "Forest_area_meta_5km.png"), forest_plot, width = 7, height = 4.5, dpi = 300, bg = "white")
print(forest_plot)

# Woody vegetation turnover by year ----
## Area of each WVSC change class within the dissolved 5 km buffer; removal and degradation are losses, establishment and enhancement gains
woody_change <- map(Wvsc_years, \(year) {
  wvsc <- read_in_area(wvsc_paths[[as.character(year)]], buffer_union)
  cell_ha <- terra::cellSize(wvsc, unit = "ha", mask = TRUE)
  terra::zonal(cell_ha, wvsc, fun = "sum", na.rm = TRUE) %>%
    set_names(c("Code", "Area_ha")) %>%
    mutate(Year = year, Buffer_ha = terra::global(cell_ha, "sum", na.rm = TRUE)[[1]])
}) %>%
  list_rbind() %>%
  filter(as.character(Code) %in% names(Wvsc_classes)) %>%
  complete(nesting(Year, Buffer_ha), Code = as.numeric(names(Wvsc_classes)), fill = list(Area_ha = 0)) %>%
  mutate(
    Class = factor(Wvsc_classes[as.character(Code)], levels = Wvsc_classes),
    Direction = if_else(Class %in% c("Removal", "Degradation"), "Loss", "Gain"),
    Area_pct = round(100 * Area_ha / Buffer_ha, 3),
    Area_ha = round(Area_ha, 1)
  ) %>%
  select(Year, Class, Direction, Area_ha, Area_pct) %>%
  arrange(Year, Class)
write_csv(woody_change, Woody_change_csv)

## Three panels, losses drawn below zero and gains above: turnover of woody area (removal, cleared to zero, vs establishment, new from zero); to its left, Hansen's annual forest loss for 2001–2012 on the same y scale, with the undated 2000–2012 gain shown as its annual mean, extending the record before WVSC begins; below, structural change in woody vegetation that persists (degradation vs enhancement), on its own y scale since growth puts many persisting pixels in enhancement every year, so it is not comparable in area to turnover
## The dashed line separates the less precise 2010–2014 WVSC maps
hansen_turnover <- tibble(Year = 2001:2012) %>%
  mutate(
    `Forest loss` = -map_dbl(Year, \(year) sum(forest_2000$Area_ha[forest_2000$Lossyear == year - 2000])),
    `Forest gain (annual mean)` = hansen_gain_ha / 12
  ) %>%
  pivot_longer(-Year, names_to = "Class", values_to = "Signed_ha")
signed_change <- woody_change %>%
  mutate(Signed_ha = if_else(Direction == "Loss", -Area_ha, Area_ha))
## Insets in panels a and b show the overall trend as cumulative net change in area since the start of each record, on one shared y range so the two are directly comparable: Hansen forest relative to 2000 with the 2000–2012 gain (straight line, as in the forest-area figure); WVSC establishment minus removal summed over years, i.e. net change in woody area (enhancement and degradation are not summed, since the same persisting pixels recur in them every year)
hansen_net <- hansen_forest %>%
  filter(Year %in% c(2000, 2012)) %>%
  mutate(Net_ha = Forest_ha - Forest_ha[1])
potapov_net <- signed_change %>%
  filter(Class %in% c("Removal", "Establishment")) %>%
  group_by(Year) %>%
  summarise(Net_ha = sum(Signed_ha)) %>%
  add_row(Year = min(Wvsc_years) - 1, Net_ha = 0, .before = 1) %>%
  mutate(Net_ha = cumsum(Net_ha))
inset_y_limits <- range(c(0, hansen_net$Net_ha, potapov_net$Net_ha))
hansen_inset <- net_change_inset(hansen_net, "Net change since 2000", sprintf("Net change = %s ha", format(round(last(hansen_net$Net_ha)), big.mark = ",")), seq(2000, 2012, 4), inset_y_limits)
potapov_inset <- net_change_inset(potapov_net, "Net change since 2009", sprintf("Establishment - removal = %s ha", format(round(last(potapov_net$Net_ha)), big.mark = ",")), seq(2010, 2024, 2), inset_y_limits)
## The shared y axis of panels a and b is extended upwards to make room for the insets above the bars
bar_range <- range(c(0, hansen_turnover$Signed_ha, filter(signed_change, Class %in% c("Removal", "Establishment"))$Signed_ha))
turnover_limits <- c(bar_range[1], bar_range[2] + 1.1 * diff(bar_range))
inset_y <- c(bar_range[2] + 0.12 * diff(bar_range), turnover_limits[2])
## Panel titles state each panel's main result; the data sources are cited in the corner note
hansen_panel <- change_panel(hansen_turnover, "No net forest loss", c("Forest loss" = "#67001f", "Forest gain (annual mean)" = "#4d9221"), seq(2002, 2012, 5), turnover_limits, seq(-500, 1000, 500)) +
  annotation_custom(ggplotGrob(hansen_inset), xmin = 2000.4, xmax = 2012.6, ymin = inset_y[1], ymax = inset_y[2]) +
  labs(x = NULL) +
  theme(legend.position = "top")
turnover_panel <- change_panel(filter(signed_change, Class %in% c("Removal", "Establishment")), "Woody cover >2 m increased", c(Removal = "#b2182b", Establishment = "#1b7837"), seq(2010, 2024, 2), turnover_limits, seq(-500, 1000, 500)) +
  annotate("segment", x = 2014.5, xend = 2014.5, y = -Inf, yend = bar_range[2], linetype = "dashed", colour = "grey40") +
  annotation_custom(ggplotGrob(potapov_inset), xmin = 2009.4, xmax = 2024.6, ymin = inset_y[1], ymax = inset_y[2]) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "top")
structure_panel <- change_panel(filter(signed_change, Class %in% c("Degradation", "Enhancement")), "Woody enhancement outpaces degradation", c(Degradation = "#f4a582", Enhancement = "#a6dba0"), seq(2010, 2024, 2)) +
  geom_vline(xintercept = 2014.5, linetype = "dashed", colour = "grey40") +
  theme(legend.direction = "horizontal")
change_plot <- plot_grid(
  ggdraw() + draw_label("Landscape context within 5 km of survey points", size = 15, x = 0.01, hjust = 0),
  plot_grid(hansen_panel, turnover_panel, nrow = 1, rel_widths = c(0.55, 1), align = "h", axis = "tb", labels = c("a", "b")),
  plot_grid(structure_panel, labels = "c"),
  ncol = 1, rel_heights = c(0.07, 1.3, 1)
) +
  draw_label("a) from Hansen et al. (2013)\nb) and c) from Potapov et al. (2025)", x = 0.99, y = 0.01, hjust = 1, vjust = 0, size = 9, colour = "grey30")
ggsave(file.path(Fig_dir, "Woody_change_meta_5km.png"), change_plot, width = 8.5, height = 8, dpi = 300, bg = "white")
print(change_plot)

# Console report ----
cat(sprintf("Meta points: %d (%d farms); dissolved %g km buffer: %.0f ha\n", nrow(meta_pts), n_distinct(meta_pts$Farm_name), Buffer_m / 1000, hansen_buffer_ha))
hansen_pct <- \(year) 100 * hansen_forest$Forest_ha[hansen_forest$Year == year] / hansen_buffer_ha
cat(sprintf("Hansen forest (with the 2000–2012 gain): %.1f%% in 2000, %.1f%% in 2012, %.1f%% in 2024; loss 2001–2012 %.0f ha, 2013–2024 %.0f ha; gain 2000–2012 %.0f ha\n",
  hansen_pct(2000), hansen_pct(2012), hansen_pct(2024), sum(forest_2000$Area_ha[forest_2000$Lossyear %in% 1:12]), sum(forest_2000$Area_ha[forest_2000$Lossyear %in% 13:24]), hansen_gain_ha))
print(woody_forest %>% mutate(Forest_pct = round(100 * Forest_ha / Buffer_ha, 1)))
