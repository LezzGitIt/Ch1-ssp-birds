## Meta sampling map ##
# Point-count sampling map zoomed to the Meta department (the sole study region for this
# chapter), on an elevation background - adapted from the full multi-ecoregion version in
# ../Ssp-bird-data-wrangling/Data_paper/Scripts/Data_paper_ecology_figs.R. That script (and this
# one's previous version) filtered to Ecoregion == "Piedemonte" rather than Departamento ==
# "Meta" - in this dataset the two happen to select the exact same 356 rows, but Departamento is
# the semantically correct field for a map that's explicitly about Meta, so this version filters
# on that instead. This repo has no Data/ folder of its own yet, so point/spatial layers are read
# directly from the sibling Ssp-bird-data-wrangling repo.

# Libraries -----------------------------------------------------------------
library(tidyverse)
library(sf)
library(cowplot)
library(ggrepel)
library(rnaturalearthdata)
library(rnaturalearth)
library(geodata)
library(terra)
library(tidyterra)
library(ggspatial)
library(conflicted)
library(patchwork)
library(magick)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)

# Round `x` to the nearest multiple of `base` - used below to bin point-count coordinates into
# grid cells for overlap counting/jittering (Themes_funs.R, sourced by the sibling script this
# was adapted from, isn't present on this machine anymore, so defined directly here instead).
mround <- function(x, base) base * round(x / base)

# Load data -------------------------------------------------------------------
wrangling_repo <- "../Ssp-bird-data-wrangling"
load(file.path(wrangling_repo, "Rdata/NE_layers_Colombia.Rdata")) # rivers_co, neCol, neColDepts, SA
Pc_locs_dc_sf <- st_read(file.path(wrangling_repo, "Derived_geospatial/shp/Pc_locs_dc.gpkg"), quiet = TRUE)

# Elevation background --------------------------------------------------------
# elevation_30s() reads from the shared cache below if already downloaded (it is, here) rather
# than re-fetching - path is a sibling of both repos, same relative depth from either one.
Envi_path <- "../Geospatial_data/Environmental"
Elev_1km <- geodata::elevation_30s(country = "Colombia", path = Envi_path)
ColElev_df <- terra::as.data.frame(Elev_1km, xy = TRUE) %>% tibble()

# Two elevation color scales sharing one raster layer: untransformed (linear) for the tightly
# zoomed Meta inset, where the whole visible area is lowland foothill and a log transform would
# over-emphasize small variation there; log-transformed for the full-Colombia context map, where
# elevation spans sea level to the Andes and the log scale is what makes the lowlands legible at
# all (otherwise the Andes alone would dominate the color range).
Col_alt_base <- ggplot() +
  geom_raster(data = ColElev_df, aes(x = x, y = y, fill = COL_elv_msk)) +
  labs(x = "Longitude", y = "Latitude", fill = "Elevation") +
  guides(fill = "none") +
  theme(axis.title = element_blank())

Col_alt_map <- Col_alt_base + scale_fill_viridis_c()
Col_alt_map_log <- Col_alt_base + scale_fill_viridis_c(trans = "log")

# Point formatting ----------------------------------------------------------
# Round to 0.03-degree cells (much finer than the 0.25 degree used for the full-country map,
# appropriate now that the view is zoomed to one department instead of all of Colombia), count
# overlapping point counts per cell, then jitter proportionally to how many share a cell so nearby
# but distinct sites don't render as one dot.
Pc_locs_meta <- Pc_locs_dc_sf %>%
  filter(Departamento == "Meta", Uniq_db != "Cipav mbd") %>%
  st_drop_geometry() %>%
  tibble() %>%
  mutate(Latitud_rd = mround(Latitud, .015), Longitud_rd = mround(Longitud, .015)) %>%
  count(Uniq_db, Ecoregion, Latitud_rd, Longitud_rd, sort = TRUE) %>%
  st_as_sf(coords = c("Longitud_rd", "Latitud_rd"), crs = 4326, remove = FALSE) %>%
  mutate(dup_count = n(), .by = c(Latitud_rd, Longitud_rd))

# Jitter amount scales with how many duplicate-count groups are actually present (rather than a
# hardcoded-length vector), so this doesn't break if Meta's data has more/fewer overlap levels
# than the 4 the original full-country version happened to have.
n_dup_groups <- n_distinct(Pc_locs_meta$dup_count)
jitter_amt <- seq(.02, .08, length.out = n_dup_groups)
Pc_locs_meta_jit <- Pc_locs_meta %>%
  group_split(dup_count) %>%
  map2(jitter_amt, \(locs, jam) st_jitter(locs, jam)) %>%
  dplyr::bind_rows()

# Inset elevation map -----------------------------------------------------
bbox_meta <- st_bbox(Pc_locs_meta_jit) + c(-.3, -.3, .3, .3)

# River names are shown on the Colombia context map below instead, not here - at this zoom every
# river label ended up either clipped at the crop edge or (once fixed) crowding the small panel;
# the line itself is enough context in the tightly zoomed inset.
meta_map <- Col_alt_map +
  geom_sf(data = Pc_locs_meta_jit, color = "red", aes(shape = Uniq_db, size = n, alpha = desc(n))) +
  geom_sf(data = rivers_co, color = "blue") +
  coord_sf(
    xlim = c(bbox_meta[1], bbox_meta[3]), ylim = c(bbox_meta[2], bbox_meta[4]),
    label_axes = "____", expand = TRUE
  ) +
  annotation_scale(location = "bl") +
  # CIPAV excluded above, so its level (which sorted first, alphabetically) is dropped from these
  # 5 remaining labels rather than the original 6
  scale_shape_discrete(
    labels = c("GAICA\ndistancia", "GAICA", "UBC & GAICA", "UBC", "Universidad de \nlos Llanos"),
    solid = FALSE
  ) +
  scale_size_continuous(range = c(3, 7)) +
  scale_alpha_continuous(range = c(.5, 1)) +
  # override.aes: points on the map itself stay red (the actual data color), but the legend key
  # swatches are drawn in black instead - deliberately different from the plotted color, per
  # request, so the legend reads clearly regardless of the red-on-viridis contrast in the map.
  guides(
    alpha = "none",
    size = guide_legend(title = "Number of \npoint counts", override.aes = list(color = "black")),
    shape = guide_legend(title = "Data collector", override.aes = list(color = "black"))
  )
meta_map

# Colombia context map ------------------------------------------------------
# The zoomed sampling map above loses all sense of where Meta sits within Colombia - this adds
# that context back: full-country department outlines (from rnaturalearth's bundled medium-
# resolution data, no download needed), Meta highlighted, and a rectangle marking exactly where
# the zoomed map above is cropped to, so the two panels are legible together, not just
# side-by-side coincidence.
co_depts <- ne_states(country = "Colombia", returnclass = "sf")
co_bbox <- st_bbox(co_depts)
meta_dept <- co_depts %>% filter(name == "Meta")

# The 4 rivers the original full-country script singled out (Magdalena, Cauca, Guainía, Meta) -
# at this whole-of-Colombia scale all 4 are real, legible, worth labeling rivers, unlike at the
# tight Meta-only zoom where Guainía turned out to be a naming collision, not the actual river.
rivers_named <- rivers_co %>% filter(name %in% c("Magdalena", "Cauca", "Guainía", "Meta"))

colombia_map <- Col_alt_map_log +
  geom_sf(data = co_depts, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf(data = meta_dept, fill = NA, color = "red", linewidth = 1) +
  geom_sf(data = rivers_co, color = "blue") +
  geom_text_repel(
    data = rivers_named,
    aes(label = name, geometry = geometry),
    stat = "sf_coordinates",
    size = 4, color = "red", max.overlaps = Inf
  ) +
  annotate(
    "rect",
    xmin = bbox_meta["xmin"], xmax = bbox_meta["xmax"],
    ymin = bbox_meta["ymin"], ymax = bbox_meta["ymax"],
    fill = NA, color = "black", linewidth = 0.6
  ) +
  coord_sf(label_axes = "____")
colombia_map

# Export ------------------------------------------------------------------
# Three separate pieces, sized/arranged manually in PowerPoint from here rather than pre-combined
# in R (connector lines between the rectangle above and the inset below are much easier to draw
# by hand once than to get right programmatically - see prior discussion). meta_map_inset drops
# the legend (it's exported on its own instead, immediately below) and is saved at the same
# width/height as colombia_map so the two start from matching scale before any manual resizing.
#
# ggsave() alone leaves a wide white margin around irregularly-shaped content like a country
# outline or a legend (whichever dimension - width or height - the content is relatively narrower
# in gets padded out to fill the requested canvas), which was requiring a manual crop in
# PowerPoint every time. save_trimmed() saves normally, then uses magick::image_trim() to strip
# any uniform white border from the saved file directly - simpler and more robust than hand-
# computing each plot's exact data aspect ratio, which wouldn't have touched the legend's own
# padding anyway.
save_trimmed <- function(filename, plot, ...) {
  ggsave(filename, plot, ...)
  magick::image_read(filename) %>%
    magick::image_trim(fuzz = 2) %>%
    magick::image_write(filename)
}

map_width <- 8
map_height <- 8
save_trimmed("Figures/Map_sampling/Map_sampling_meta_colombia.png", colombia_map, bg = "white", width = map_width, height = map_height)

meta_map_inset <- meta_map + theme(legend.position = "none")
save_trimmed("Figures/Map_sampling/Map_sampling_meta_inset.png", meta_map_inset, bg = "white", width = map_width, height = map_height)

# cowplot::get_legend() pulls the legend out as its own grob/gtable (not a ggplot object, so
# plain ggplot()/ggsave() on it directly won't work) - cowplot::ggdraw() wraps it back into a
# ggplot-drawable object so ggsave() can render it to its own file.
meta_map_legend <- cowplot::get_legend(meta_map)
save_trimmed("Figures/Map_sampling/Map_sampling_meta_legend.png", cowplot::ggdraw(meta_map_legend), bg = "white", width = 4, height = 5)
