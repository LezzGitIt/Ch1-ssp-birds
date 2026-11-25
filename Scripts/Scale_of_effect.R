## PhD birds in silvopastoral landscapes ##
# Scale of effect of landscape canopy cover on occupancy: extracts woody canopy cover (trees / shrubs > 2 m, Colombia Woody Vegetation Structure and Change, ~27 m) in rings around each point count over a wide range of radii, matched to the survey year, then checks how the covariate behaves across radii (variation, correlation among radii, collinearity with the other occupancy covariates) before any model is fitted.

### Rings, not discs: each ring excludes the inner Ring_inner_m around the point, so the landscape measure does not duplicate the point's own Habitat (as in Karp et al. 2019).
### Only the surveyed location x period rows of the stacked data list (the cells that inform the multi-season model) enter the diagnostics.
### Scale-of-effect fits (last section, run only when named on the command line): the primary model with one ring's canopy cover added, one fit per radius, plus the primary model refitted unchanged as the baseline. The fits made 2026-11-07 (baseline, r150, r300, r500) used the earlier primary model (signed distance, Year1, loadings started at 0); they were superseded and deleted 2026-11-09.
### The habitat-by-landscape-cover and habitat-conversion models built on the 1.5 km ring were fitted here (as soe_<model>, 2026-11-21 to 23) and moved to the multi-season model script on 2026-11-24, which reads this script's per-row covariate table.

# Load libraries ----
library(tidyverse)
library(cowplot)
ggplot2::theme_set(cowplot::theme_cowplot())
source("Scripts/spOcc_fns.R")

# Configuration ----
## Outer ring radii (m), roughly even on a log scale from about one territory to well beyond the ~1 km median scale of effect for birds (Jackson & Fahrig 2015); 10 km is the widest the Farm_diversity raster clips cover
## 100 and 200 m dropped (user, 2026-11-07): ring cover there is nearly redundant with neighbouring radii and largely duplicates Habitat and Dist_forest (r with Dist_forest about -0.75, Habitat R² up to 0.59)
Radii_m <- c(150, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000)
Ring_inner_m <- 50

## Inputs: the stacked list the multi-season list was built from (one row per location x year-group x season), and the wrangling repo's survey and site tables
msOcc_l_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm.rds"
Wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"

## Canopy rasters: farm x year clips made by Farm_diversity's canopy-buffer script (10.5 km around each farm's points); a missing clip is cropped from the national raster instead
Clip_dir <- "../Farm_diversity/Data/Geospatial/canopy_rasters"
National_dir <- "../Geospatial_data/Environmental"
Raster_years <- c(2013, 2014, 2016, 2017, 2019, 2022, 2024)

## Rings wider than the farm clips reach (user, 2026-11-15: migrants may perceive the landscape at large scales, but the study extent, at most 66 km between points, leaves little contrast beyond this), from the national canopy maps into their own table, so the farm-scale ring table stays as the 98-species fits used it
Wide_radii_m <- 20000

## Canopy height (user, 2026-11-17: cover alone cannot separate low woody vegetation such as live fences, scrub and silvopasture trees from tall forest): the height maps (Colombia_WVCH, decimetres, 0 exactly where cover is 0) give the mean height of the woody canopy in each ring and the share of each ring under canopy at or above each threshold (m)
Tall_canopy_m <- c(10, 15)
## Radii of the cover-by-woody-canopy-height fits ('h<radius>', user 2026-11-17): from 500 m, as silvopasture cover within 300 m (ssp) and the distance terms already describe the immediate surroundings, to 2.5 km, beyond which woody canopy height correlates ≥ 0.79 with ring cover; 3 and 5 km added (user, 2026-11-20) to find where the interaction can no longer be separated from cover (r 0.79 and 0.85)
Height_fit_radii_m <- c(500, 1000, 1500, 2000, 2500, 3000, 5000)

## Outputs
Rings_csv <- "Derived/Excels/Canopy_rings_by_location.csv"
Rings_wide_csv <- "Derived/Excels/Canopy_rings_by_location_20km.csv"
Height_rings_csv <- "Derived/Excels/Canopy_height_rings_by_location.csv"
Fit_covs_csv <- "Derived/Excels/Scale_of_effect_fit_covs.csv"
Fig_dir <- "Figures/Scale_of_effect"

## Fits to run, from the command line: 'baseline', 'r<radius>' (e.g. r300), 'h<radius>' / 'ha<radius>' (cover-by-height fits), or 'smoke' (10 batches, 1 chain, with the 1 km ring; checks the model terms); none when the script is sourced
run_sections <- commandArgs(trailingOnly = TRUE)


## Primary model (Multi_spp_100 row 107) with the mean-survey-date coding in place of Year1: distance to the forest edge hinge-coded (Dist_open outside forest, Dist_inside inside it), 2 spatial factors with the loadings started at row 90's posterior means; every fit here is this model plus at most one ring term, all with the same seed
tmsOcc_l_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_98_hab_fix.rds"
Lambda_inits_path <- "Derived/Rdata/Model_output/Temporal/lambda_inits_prec_nofarm_long_f2.rds"
Primary_occ_formula <- ~ Habitat * Tot_prec + Dist_open + Dist_inside + ssp + Year_cont
## MCMC: 24,000 iterations with 12,000 burn-in, thinned by 30 (400 kept per chain, the most the species x site x period arrays allow in memory), 40% of the 60,000 iterations of row 107; row 107 refitted at this length matched it (all community coefficients within 0.08 posterior SD, max community Rhat 1.04, minimum community ESS 198; 2026-11-09)
Primary_fit_settings <- list(
  det.formula = ~ Forest_bin + Pc_start + I(Pc_start^2) + Prob_breed + Sampling_day_log + (1 | Collector_team_num),
  priors = list(phi.unif = list(3 / 25e3, 3 / 23e3)),
  n.neighbors = 5, n.factors = 2, batch.length = 25, ar1 = FALSE,
  mcmc = list(n.batch = 960, n.burn = 12000, n.thin = 30, n.chains = 3),
  out_dir = "Derived/Rdata/Model_output/Temporal/Scale_of_effect",
  seed = 20261107
)
## Cover-by-height refits at full length ('ha<radius>', user 2026-11-20 night): habitat and rainfall additive, as the h fits should have had them, and the 60,000 iterations of row 107 (30,000 burn-in, thinned by 75, still 400 kept per chain), since the 1.5 km h fit reached a max community Rhat of 1.14 at 24,000
Additive_occ_formula <- update(Primary_occ_formula, ~ . - Habitat:Tot_prec)
Mcmc_long <- list(n.batch = 2400, n.burn = 30000, n.thin = 75, n.chains = 3)

# Functions ----
# Ring means from exactextractr's coverage-weighted disc sums and cell counts (one row per location x raster year x radius, the inner disc included): ring mean = (outer sum - inner sum) / (outer count - inner count)
ring_means <- function(disc_sums, inner_m) {
  disc_sums %>%
    filter(Radius_m != inner_m) %>%
    left_join(disc_sums %>% filter(Radius_m == inner_m) %>% select(Id_survey_no_dc, Raster_year, sum_inner = sum, count_inner = count), by = c("Id_survey_no_dc", "Raster_year")) %>%
    mutate(Canopy_cover_ring = round((sum - sum_inner) / (count - count_inner), 3)) %>%
    select(Id_survey_no_dc, Raster_year, Radius_m, Canopy_cover_ring) %>%
    arrange(Id_survey_no_dc, Raster_year, Radius_m)
}

# Load data ----
## Surveyed rows of the model, with the survey year taken as in the wrangling script (latest year within a location x year-group x season) and snapped to the nearest raster year (2025-26 surveys use 2024)
Event_covs <- read_csv(paste0(Wrangling_excels, "Event_covs.csv"), show_col_types = FALSE)
Event_years <- Event_covs %>%
  distinct(Id_survey, Id_survey_no_dc, Year_grp, Season, Year) %>%
  slice_max(by = c(Id_survey, Year_grp, Season), order_by = Year, with_ties = FALSE)

Site_locs <- read_csv(paste0(Wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  select(Id_survey_no_dc, Id_scr, Long, Lat)

## The data list stores Id_survey_no_dc and Id_scr as numeric codes, so the location and farm come from Id_survey via the wrangling tables; the location code is kept as Site_code, the site name in the multi-season list
Rows <- readRDS(msOcc_l_path)$occ.covs %>%
  select(Id_survey, Year_grp, Season, Site_code = Id_survey_no_dc, Habitat, Dist_forest, Canopy_cover) %>%
  mutate(across(c(Id_survey, Year_grp, Season, Site_code), as.character)) %>%
  left_join(Event_years, by = c("Id_survey", "Year_grp", "Season")) %>%
  left_join(Site_locs, by = "Id_survey_no_dc") %>%
  mutate(Raster_year = map_dbl(Year, \(year) Raster_years[which.min(abs(Raster_years - year))]))
stopifnot(!anyNA(Rows$Year), !anyNA(Rows$Long))

Loc_years <- Rows %>% distinct(Id_survey_no_dc, Id_scr, Long, Lat, Raster_year)

# Extract canopy cover in rings ----
### Run once, then read from the CSV. Ring mean = (canopy sum in the outer disc - sum in the inner disc) / (cell count outer - cell count inner), with exactextractr's coverage-weighted sums and counts, so every ring is exactly its own disc minus that point's inner disc (geometric erasing would also cut out neighbouring points' inner discs).
if (file.exists(Rings_csv)) {
  Canopy_rings <- read_csv(Rings_csv, show_col_types = FALSE)
  Canopy_rings <- filter(Canopy_rings, Radius_m %in% Radii_m)
  stopifnot(setequal(unique(Canopy_rings$Radius_m), Radii_m), nrow(anti_join(Loc_years, Canopy_rings, by = c("Id_survey_no_dc", "Raster_year"))) == 0)
} else {
  Disc_sums <- Loc_years %>%
    group_split(Id_scr, Raster_year) %>%
    map(\(grp) {
      pts <- terra::vect(grp, geom = c("Long", "Lat"), crs = "EPSG:4326")
      clip_path <- file.path(Clip_dir, sprintf("farm_%s_%d.tif", grp$Id_scr[1], grp$Raster_year[1]))
      canopy_r <- if (file.exists(clip_path)) {
        terra::rast(clip_path)
      } else {
        message("No Farm_diversity clip for farm ", grp$Id_scr[1], " / ", grp$Raster_year[1], "; cropping the national raster")
        terra::crop(terra::rast(file.path(National_dir, sprintf("WVCC/Colombia_WVCC_%d_50km.tif", grp$Raster_year[1]))), terra::buffer(pts, max(Radii_m) + 500), snap = "out")
      }
      # The widest disc must lie inside the raster
      need <- as.vector(terra::ext(terra::buffer(pts, max(Radii_m))))
      have <- as.vector(terra::ext(canopy_r))
      stopifnot(need[c("xmin", "ymin")] >= have[c("xmin", "ymin")], need[c("xmax", "ymax")] <= have[c("xmax", "ymax")])
      map(c(Ring_inner_m, Radii_m), \(radius) {
        discs <- sf::st_as_sf(terra::buffer(pts, width = radius))
        exactextractr::exact_extract(canopy_r, discs, c("sum", "count"), progress = FALSE) %>%
          mutate(Id_survey_no_dc = grp$Id_survey_no_dc, Raster_year = grp$Raster_year, Radius_m = radius)
      }) %>% list_rbind()
    }) %>%
    list_rbind()

  Canopy_rings <- ring_means(Disc_sums, Ring_inner_m)
  write_csv(Canopy_rings, Rings_csv)
}
stopifnot(!anyNA(Canopy_rings$Canopy_cover_ring))

## Wider rings, run once: the same ring means from each raster year's national map, cropped around that year's locations
if (!file.exists(Rings_wide_csv)) {
  Loc_years %>%
    group_split(Raster_year) %>%
    map(\(grp) {
      pts <- terra::vect(grp, geom = c("Long", "Lat"), crs = "EPSG:4326")
      canopy_r <- terra::crop(terra::rast(file.path(National_dir, sprintf("WVCC/Colombia_WVCC_%d_50km.tif", grp$Raster_year[1]))), terra::buffer(pts, max(Wide_radii_m) + 500), snap = "out")
      map(c(Ring_inner_m, Wide_radii_m), \(radius) {
        exactextractr::exact_extract(canopy_r, sf::st_as_sf(terra::buffer(pts, width = radius)), c("sum", "count"), progress = FALSE) %>%
          mutate(Id_survey_no_dc = grp$Id_survey_no_dc, Raster_year = grp$Raster_year, Radius_m = radius)
      }) %>% list_rbind()
    }) %>%
    list_rbind() %>%
    ring_means(Ring_inner_m) %>%
    write_csv(Rings_wide_csv)
}

# Extract canopy height in rings ----
### Run once, for every radius including the widest, from each raster year's national cover and height maps cropped around that year's locations. Per pixel, cover share c (0-1) and height h (m); per disc, coverage-weighted sums of c, c x h and c x [h >= threshold], differenced to rings as for cover. Woody canopy height = ring sum of c x h / ring sum of c (mean height weighted by each pixel's cover; undefined in a ring with no canopy); tall-canopy cover = ring sum of c x [h >= threshold] / ring cell count; ring cover is recomputed the same way as a check against the cover rings.
if (file.exists(Height_rings_csv)) {
  Height_rings <- read_csv(Height_rings_csv, show_col_types = FALSE)
  stopifnot(all(c(Radii_m, Height_fit_radii_m, Wide_radii_m) %in% Height_rings$Radius_m), nrow(anti_join(Loc_years, Height_rings, by = c("Id_survey_no_dc", "Raster_year"))) == 0)
} else {
  Height_rings <- Loc_years %>%
    group_split(Raster_year) %>%
    map(\(grp) {
      pts <- terra::vect(grp, geom = c("Long", "Lat"), crs = "EPSG:4326")
      window <- terra::buffer(pts, max(c(Radii_m, Wide_radii_m)) + 500)
      cover_r <- terra::crop(terra::rast(file.path(National_dir, sprintf("WVCC/Colombia_WVCC_%d_50km.tif", grp$Raster_year[1]))), window, snap = "out") / 100
      height_r <- terra::crop(terra::rast(file.path(National_dir, sprintf("WVCH/Colombia_WVCH_%d_50km.tif", grp$Raster_year[1]))), window, snap = "out") / 10
      layers <- terra::rast(c(list(Cover = cover_r, Cover_height = cover_r * height_r), map(Tall_canopy_m, \(m) cover_r * (height_r >= m)) %>% set_names(paste0("Tall_cover_", Tall_canopy_m, "m"))))
      map(sort(unique(c(Ring_inner_m, Radii_m, Height_fit_radii_m, Wide_radii_m))), \(radius) {
        exactextractr::exact_extract(layers, sf::st_as_sf(terra::buffer(pts, width = radius)), c("sum", "count"), progress = FALSE) %>%
          select(starts_with("sum."), count = count.Cover) %>%
          mutate(Id_survey_no_dc = grp$Id_survey_no_dc, Raster_year = grp$Raster_year, Radius_m = radius)
      }) %>% list_rbind()
    }) %>%
    list_rbind() %>%
    pivot_longer(c(starts_with("sum."), count), names_to = "Stat") %>%
    mutate(value = value - value[Radius_m == Ring_inner_m], .by = c(Id_survey_no_dc, Raster_year, Stat)) %>%
    filter(Radius_m != Ring_inner_m) %>%
    pivot_wider(names_from = Stat) %>%
    mutate(
      Canopy_cover_check = sum.Cover / count,
      Canopy_height_woody_ring = if_else(sum.Cover > 0, sum.Cover_height / sum.Cover, NA_real_),
      across(starts_with("sum.Tall_cover_"), \(x) x / count)
    ) %>%
    rename_with(\(x) paste0(sub("sum.", "", x, fixed = TRUE), "_ring"), starts_with("sum.Tall_cover_")) %>%
    select(Id_survey_no_dc, Raster_year, Radius_m, Canopy_cover_check, Canopy_height_woody_ring, starts_with("Tall_cover_")) %>%
    mutate(across(-c(Id_survey_no_dc, Raster_year, Radius_m), \(x) round(x, 3))) %>%
    arrange(Id_survey_no_dc, Raster_year, Radius_m)
  write_csv(Height_rings, Height_rings_csv)
}

# Covariate diagnostics across radii ----
Rows_rings <- Rows %>%
  left_join(Canopy_rings, by = c("Id_survey_no_dc", "Raster_year"), relationship = "many-to-many")

## Per radius: spread (SD, coefficient of variation = SD / mean), share of the variance lying among farms (R² of a farm-only linear model), correlation with Dist_forest, variance explained by Habitat, and correlation with the model's current 1 km Canopy_cover disc as a check on the extraction
Scale_summary <- Rows_rings %>%
  summarize(
    .by = Radius_m,
    Mean = mean(Canopy_cover_ring),
    SD = sd(Canopy_cover_ring),
    CV = SD / Mean,
    Among_farm_share = summary(lm(Canopy_cover_ring ~ factor(Id_scr)))$r.squared,
    r_Dist_forest = cor(Canopy_cover_ring, Dist_forest),
    R2_Habitat = summary(lm(Canopy_cover_ring ~ Habitat))$r.squared,
    r_Canopy_cover_1km_disc = cor(Canopy_cover_ring, Canopy_cover)
  ) %>%
  arrange(Radius_m)

Radius_cor <- Rows_rings %>%
  select(Id_survey, Year_grp, Season, Radius_m, Canopy_cover_ring) %>%
  pivot_wider(names_from = Radius_m, values_from = Canopy_cover_ring) %>%
  select(as.character(Radii_m)) %>%
  cor()

## Canopy height variables per radius: correlation with ring cover (whether height adds anything beyond the amount of canopy), with Dist_forest, and variance explained by Habitat; rings without canopy are dropped for woody canopy height and counted
Rows_height <- Rows %>%
  left_join(Height_rings, by = c("Id_survey_no_dc", "Raster_year"), relationship = "many-to-many")
Cover_check <- inner_join(Rows_height, Rows_rings, by = c("Id_survey", "Year_grp", "Season", "Radius_m"))
stopifnot(cor(Cover_check$Canopy_cover_check, Cover_check$Canopy_cover_ring) > 0.999)

Height_summary <- Rows_height %>%
  pivot_longer(c(Canopy_height_woody_ring, starts_with("Tall_cover_")), names_to = "Variable") %>%
  filter(!is.na(value)) %>%
  summarize(
    .by = c(Variable, Radius_m),
    N_rows = n(),
    Mean = mean(value),
    SD = sd(value),
    r_Canopy_cover_ring = cor(value, Canopy_cover_check),
    r_Dist_forest = cor(value, Dist_forest),
    R2_Habitat = summary(lm(value ~ Habitat))$r.squared
  ) %>%
  arrange(Variable, Radius_m)

# Figures ----
dir.create(Fig_dir, recursive = TRUE, showWarnings = FALSE)
Radius_axis <- scale_x_log10(breaks = Radii_m, labels = \(x) ifelse(x >= 1000, paste0(x / 1000, " km"), paste0(x, " m")))

## Variation in canopy cover across radii: larger rings average over more of the region, shrinking the spread that lets a scale be detected (Miguet et al. 2016)
variation_plot <- Scale_summary %>%
  select(Radius_m, `Coefficient of variation` = CV, `Share of variance among farms` = Among_farm_share) %>%
  pivot_longer(-Radius_m) %>%
  ggplot(aes(Radius_m, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  Radius_axis +
  labs(x = paste0("Ring outer radius (inner ", Ring_inner_m, " m excluded)"), y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(Fig_dir, "Canopy_rings_variation.png"), variation_plot, width = 7, height = 7, bg = "white")
print(variation_plot)

## Correlation of canopy cover between radii
radius_cor_plot <- as_tibble(Radius_cor, rownames = "Radius_a") %>%
  pivot_longer(-Radius_a, names_to = "Radius_b", values_to = "r") %>%
  mutate(across(c(Radius_a, Radius_b), \(x) factor(x, levels = as.character(Radii_m)))) %>%
  ggplot(aes(Radius_a, Radius_b, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.6) +
  scale_fill_viridis_c(limits = c(min(0, min(Radius_cor)), 1)) +
  labs(x = "Ring outer radius (m)", y = "Ring outer radius (m)", fill = "Pearson r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(Fig_dir, "Canopy_rings_correlation.png"), radius_cor_plot, width = 8, height = 7, bg = "white")
print(radius_cor_plot)

## Collinearity with the other occupancy covariates of the primary model
collinearity_plot <- Scale_summary %>%
  select(Radius_m, `Correlation with Dist_forest (r)` = r_Dist_forest, `Variance explained by Habitat (R²)` = R2_Habitat) %>%
  pivot_longer(-Radius_m) %>%
  ggplot(aes(Radius_m, value)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line() +
  geom_point() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  Radius_axis +
  labs(x = paste0("Ring outer radius (inner ", Ring_inner_m, " m excluded)"), y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(Fig_dir, "Canopy_rings_collinearity.png"), collinearity_plot, width = 7, height = 7, bg = "white")
print(collinearity_plot)

## Canopy height variables against ring cover and the other occupancy covariates
height_plot <- Height_summary %>%
  select(Variable, Radius_m, `Correlation with ring canopy cover (r)` = r_Canopy_cover_ring, `Correlation with Dist_forest (r)` = r_Dist_forest, `Variance explained by Habitat (R²)` = R2_Habitat) %>%
  pivot_longer(-c(Variable, Radius_m)) %>%
  mutate(Variable = str_remove(Variable, "_ring$")) %>%
  ggplot(aes(Radius_m, value, colour = Variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line() +
  geom_point() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = sort(unique(Height_summary$Radius_m)), labels = \(x) ifelse(x >= 1000, paste0(x / 1000, " km"), paste0(x, " m"))) +
  labs(x = paste0("Ring outer radius (inner ", Ring_inner_m, " m excluded)"), y = NULL, colour = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path(Fig_dir, "Canopy_height_rings_collinearity.png"), height_plot, width = 7, height = 9, bg = "white")
print(height_plot)

# Report ----
cat(nrow(Rows), "surveyed rows,", nrow(Loc_years), "location x raster-year combinations,", n_distinct(Rows$Id_scr), "farms\n")
print(Scale_summary %>% mutate(across(-Radius_m, \(x) round(x, 3))), width = Inf)
cat("Surveyed rows with no canopy in the ring (woody canopy height undefined), by radius:\n")
print(Rows_height %>% summarize(.by = Radius_m, N_no_canopy = sum(is.na(Canopy_height_woody_ring))) %>% filter(N_no_canopy > 0))
print(Height_summary %>% mutate(across(-c(Variable, Radius_m, N_rows), \(x) round(x, 3))), n = Inf, width = Inf)

# Fit covariates per surveyed row ----
## Mean survey date of each location-period (the stacked and multi-season lists predate the date coding and carry Year1 only)
Rows_dates <- Rows %>%
  left_join(survey_dates(mutate(Event_covs, Id_survey = as.character(Id_survey))), by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one")
stopifnot(!anyNA(Rows_dates$Survey_date))
## Original-unit ring cover (%), woody canopy height (m) and mean survey date of every surveyed row at the cover-by-height radii, for scripts that predict from these fits (e.g. marginal effects); the fits standardise each over these rows
Rows_height %>%
  filter(Radius_m %in% Height_fit_radii_m) %>%
  left_join(select(Rows_dates, Id_survey, Year_grp, Season, Survey_date), by = c("Id_survey", "Year_grp", "Season"), relationship = "many-to-one") %>%
  mutate(Cover_pct = 100 * Canopy_cover_check, Height_m = Canopy_height_woody_ring) %>%
  select(Id_survey, Year_grp, Season, Radius_m, Cover_pct, Height_m, Survey_date) %>%
  mutate(across(c(Cover_pct, Height_m), \(x) round(x, 3)), Survey_date = round(Survey_date, 5)) %>%
  arrange(Radius_m, Id_survey, Year_grp, Season) %>%
  write_csv(Fit_covs_csv)

# Scale-of-effect fits ----
### Each fit is the primary model plus one ring's canopy cover as an additive term ('r<radius>'), or plus one ring's canopy cover, woody canopy height and their interaction ('h<radius>', user 2026-11-17; both from the canopy height extraction), or the primary model unchanged ('baseline'); radii are then compared by WAIC relative to the baseline and by the canopy coefficient. No repeat run of the baseline exists to gauge how much WAIC moves between runs; the closest pair is row 107 and its shorter refit (row 113), 4 WAIC apart, but the jagged profile across neighbouring radii (2026-11-10) suggests differences of ~20 are within run-to-run noise.
### Run headless, two processes at a time (a pair of fits took about as long each as one alone, 2026-11-05), started 15-30 min apart so the memory-heavy assembly of the posterior arrays at the end does not coincide: caffeinate -is Rscript Scripts/Scale_of_effect.R baseline r150 r300 > <log> 2>&1
Soe_fit_tags <- intersect(run_sections, c("smoke", "baseline", paste0("r", Radii_m), paste0("h", Height_fit_radii_m), paste0("ha", Height_fit_radii_m)))
if (length(Soe_fit_tags) > 0) {
  ## Ring canopy cover as one site x period matrix per radius, in the multi-season list's site and period order: standardised over the surveyed rows (as the wrangling script scales every covariate), unsurveyed cells filled from the location's nearest surveyed period (as the list's other time-varying covariates)
  tmsOcc_l <- readRDS(tmsOcc_l_path)
  Site_levels <- dimnames(tmsOcc_l$y)[[2]]
  Time_levels <- dimnames(tmsOcc_l$y)[[3]]
  stopifnot(nrow(distinct(Rows, Site_code, Year_grp, Season)) == sum(apply(!is.na(tmsOcc_l$y[1, , , , drop = FALSE]), c(2, 3), any)))
  Ring_covs <- Rows_rings %>%
    mutate(Canopy_ring = as.numeric(scale(Canopy_cover_ring)), .by = Radius_m) %>%
    split(.$Radius_m) %>%
    map(\(d) site_period_matrix(d$Canopy_ring, d$Site_code, paste(d$Year_grp, d$Season, sep = "."), Site_levels, Time_levels)) %>%
    set_names(\(radius) paste0("Canopy_r", radius))
  ## Ring cover and woody canopy height for the 'h' fits, standardised the same way; a surveyed ring with no canopy (none from 500 m out) would have no woody canopy height and is set to that radius's mean (0 once standardised), leaving its zero cover to the cover term
  Height_covs <- Rows_height %>%
    filter(Radius_m %in% Height_fit_radii_m) %>%
    mutate(across(c(Canopy_cover_check, Canopy_height_woody_ring), \(x) replace_na(as.numeric(scale(x)), 0)), .by = Radius_m) %>%
    split(.$Radius_m) %>%
    imap(\(d, radius) list(d$Canopy_cover_check, d$Canopy_height_woody_ring) %>%
      map(\(x) site_period_matrix(x, d$Site_code, paste(d$Year_grp, d$Season, sep = "."), Site_levels, Time_levels)) %>%
      set_names(paste0(c("Cover_h", "Height_h"), radius))) %>%
    list_flatten(name_spec = "{inner}")
  ## Hinge distance terms, from Dist_forest back in metres (Habitat is one value per location, so it is recycled down the period columns)
  oc <- tmsOcc_l$occ.covs
  Dist_scaling <- read_csv(str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv"), show_col_types = FALSE) %>%
    filter(Process == "occupancy", Covariate == "Dist_forest")
  Hinge_covs <- hinge_distance(oc$Dist_forest * Dist_scaling$SD + Dist_scaling$Mean, matrix(oc$Habitat == "Bosque", length(Site_levels), length(Time_levels)), Dist_scaling$SD)
  ## Year_cont, the long-term time trend: mean survey date (decimal years) of each location-period (this list predates the date coding and carries Year1 only), standardised over the surveyed rows, unsurveyed cells at the period mean
  Year_cont <- site_period_matrix(as.numeric(scale(Rows_dates$Survey_date)), Rows_dates$Site_code, paste(Rows_dates$Year_grp, Rows_dates$Season, sep = "."), Site_levels, Time_levels, fill = "period_mean")
  tmsOcc_l$occ.covs <- c(oc, Ring_covs, Height_covs, Hinge_covs, list(Year_cont = Year_cont))
  Lambda_inits <- readRDS(Lambda_inits_path)
  stopifnot(identical(rownames(Lambda_inits), dimnames(tmsOcc_l$y)[[1]]))

  walk(Soe_fit_tags, \(fit_tag) {
    additive <- str_starts(fit_tag, "ha")
    base_formula <- if (additive) Additive_occ_formula else Primary_occ_formula
    ring_term <- if (fit_tag == "baseline") NULL else if (fit_tag == "smoke") "Canopy_r1000" else if (str_starts(fit_tag, "h")) paste0(c("Cover_h", "Height_h"), str_remove(fit_tag, "^ha?")) else paste0("Canopy_", fit_tag)
    occ_formula <- if (is.null(ring_term)) base_formula else update(base_formula, as.formula(paste("~ . +", paste(ring_term, collapse = " * "))))
    if (length(ring_term) == 2) ring_term <- c(ring_term, paste(ring_term, collapse = ":"))
    mod <- fit_temporal(tmsOcc_l, occ_formula, paste0("soe_", fit_tag), settings = Primary_fit_settings, mcmc = if (additive) Mcmc_long else Primary_fit_settings$mcmc, inits = list(lambda = list(Lambda_inits)), smoke = fit_tag == "smoke")
    stopifnot(setequal(colnames(mod$beta.comm.samples), c(colnames(model.matrix(base_formula, as.data.frame(map(tmsOcc_l$occ.covs, \(x) if (is.matrix(x)) x[, 1] else x)))), ring_term)))
    if (fit_tag == "smoke") cat("Smoke fit OK:", deparse(occ_formula), "\n")
  })
}

# Compare radii ----
### Reads the saved fit summaries only (not the ~1.9 GB fits), so it runs whenever the script is sourced and the baseline plus at least one radius have been fitted. ΔWAIC < 0 means the ring term improves on the baseline; the canopy coefficient is on the standardised scale (per 1 SD of that radius) and per 10 percentage points of canopy cover (that radius's SD back-transformed), so radii with different spreads can be compared.
Soe_summary_files <- list.files(Primary_fit_settings$out_dir, pattern = "^Occ_mod_98spp_soe_(baseline|r[0-9]+)_.*_summary\\.rds$", full.names = TRUE)
if (any(str_detect(Soe_summary_files, "_soe_baseline_")) && length(Soe_summary_files) > 1) {
  Soe_summaries <- map(Soe_summary_files, readRDS) %>%
    set_names(str_match(basename(Soe_summary_files), "_soe_(baseline|r[0-9]+)_")[, 2])
  stopifnot(!anyDuplicated(names(Soe_summaries)))
  Baseline_waic <- Soe_summaries$baseline$waic[["WAIC"]]

  Radius_comparison <- imap(Soe_summaries, \(s, fit_tag) {
    ring_term <- if (fit_tag == "baseline") NA_character_ else paste0("Canopy_", fit_tag)
    canopy_q <- if (is.na(ring_term)) c("50%" = NA_real_, "2.5%" = NA_real_, "97.5%" = NA_real_) else s$quantiles$beta.comm[, ring_term]
    tibble(
      Fit = fit_tag,
      Radius_m = if (is.na(ring_term)) NA_real_ else as.numeric(str_remove(fit_tag, "^r")),
      WAIC = s$waic[["WAIC"]],
      pD = s$waic[["pD"]],
      dWAIC = WAIC - Baseline_waic,
      Canopy_median = canopy_q[["50%"]],
      Canopy_lo = canopy_q[["2.5%"]],
      Canopy_hi = canopy_q[["97.5%"]],
      Canopy_P_pos = if (is.na(ring_term)) NA_real_ else s$prob_positive$beta.comm[[ring_term]],
      Canopy_tau_sq = if (is.na(ring_term)) NA_real_ else s$quantiles$tau.sq.beta["50%", ring_term],
      HabitatSsp = s$quantiles$beta.comm["50%", "HabitatSsp"],
      HabitatBosque = s$quantiles$beta.comm["50%", "HabitatBosque"],
      Dist_inside = s$quantiles$beta.comm["50%", "Dist_inside"],
      Max_comm_rhat = s$max_comm_rhat,
      Min_ess = s$min_ess[[1]],
      Loadings_median_ess = s$loadings_median_ess
    )
  }) %>%
    list_rbind() %>%
    left_join(select(Scale_summary, Radius_m, Ring_SD = SD), by = "Radius_m") %>%
    mutate(across(c(Canopy_median, Canopy_lo, Canopy_hi), \(x) x / Ring_SD * 10, .names = "{.col}_per10pp")) %>%
    arrange(!is.na(Radius_m), Radius_m) %>%
    mutate(across(where(is.double), \(x) round(x, 3)))
  write_csv(Radius_comparison, "Derived/Excels/Scale_of_effect_radius_comparison.csv")

  ## ΔWAIC and the community canopy coefficient (per 10 percentage points, 95% CrI) against ring radius
  radius_fit_plot <- Radius_comparison %>%
    filter(!is.na(Radius_m)) %>%
    select(Radius_m, dWAIC, Canopy_median_per10pp, Canopy_lo_per10pp, Canopy_hi_per10pp) %>%
    pivot_longer(c(dWAIC, Canopy_median_per10pp), names_to = "name", values_to = "value") %>%
    mutate(
      lo = if_else(name == "dWAIC", NA_real_, Canopy_lo_per10pp),
      hi = if_else(name == "dWAIC", NA_real_, Canopy_hi_per10pp),
      name = recode(name, dWAIC = "WAIC minus baseline WAIC (lower is better)", Canopy_median_per10pp = "Community canopy coefficient\n(logit, per 10 percentage points)")
    ) %>%
    ggplot(aes(Radius_m, value)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_linerange(aes(ymin = lo, ymax = hi), na.rm = TRUE) +
    geom_line() +
    geom_point() +
    facet_wrap(~name, ncol = 1, scales = "free_y") +
    Radius_axis +
    labs(x = paste0("Ring outer radius (inner ", Ring_inner_m, " m excluded)"), y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(Fig_dir, "Scale_of_effect_radii.png"), radius_fit_plot, width = 7, height = 7, bg = "white")
  print(radius_fit_plot)

  print(select(Radius_comparison, Fit, WAIC, dWAIC, pD, starts_with("Canopy"), -Canopy_tau_sq, Max_comm_rhat), width = Inf)
}

## Cover-by-height fits ----
### The 'h' fits (ring cover, woody canopy height and their interaction) against the same baseline, from their saved summaries; coefficients are community-level means on the standardised scale (per 1 SD of that radius's cover or height), with the cover-only fit at the same radius for comparison
Height_summary_files <- list.files(Primary_fit_settings$out_dir, pattern = "^Occ_mod_98spp_soe_h[0-9]+_.*_summary\\.rds$", full.names = TRUE)
if (exists("Baseline_waic") && length(Height_summary_files) > 0) {
  Height_summaries <- map(Height_summary_files, readRDS) %>%
    set_names(str_match(basename(Height_summary_files), "_soe_h([0-9]+)_")[, 2])
  Height_terms <- c(Cover = "Cover_h", Height = "Height_h", Interaction = "Interaction")
  radius_label <- \(x) ifelse(x >= 1000, paste0(x / 1000, " km"), paste0(x, " m"))
  term_name <- \(term, radius) if (term == "Interaction") paste0("Cover_h", radius, ":Height_h", radius) else paste0(term, radius)

  Height_comparison <- imap(Height_summaries, \(s, radius) map(names(Height_terms), \(term) {
    k <- term_name(Height_terms[[term]], radius)
    tibble(Radius_m = as.numeric(radius), Term = term, Median = s$quantiles$beta.comm["50%", k], Lo = s$quantiles$beta.comm["2.5%", k], Hi = s$quantiles$beta.comm["97.5%", k], P_pos = s$prob_positive$beta.comm[[k]])
  }) %>%
    list_rbind() %>%
    mutate(dWAIC = s$waic[["WAIC"]] - Baseline_waic, Max_comm_rhat = s$max_comm_rhat, Iterations = s$mcmc[[1]])) %>%
    list_rbind() %>%
    left_join(select(Radius_comparison, Radius_m, dWAIC_cover_only = dWAIC), by = "Radius_m") %>%
    arrange(Radius_m, factor(Term, names(Height_terms))) %>%
    mutate(across(where(is.double), \(x) round(x, 3)))
  write_csv(Height_comparison, "Derived/Excels/Scale_of_effect_height_comparison.csv")

  ### ΔWAIC of both fit types and the three community coefficients against ring radius
  height_radii_plot <- bind_rows(
    distinct(Height_comparison, Radius_m, `Cover x height` = dWAIC, `Cover only` = dWAIC_cover_only) %>%
      pivot_longer(-Radius_m, names_to = "Fit", values_to = "Median") %>%
      mutate(Panel = "WAIC minus baseline WAIC (lower is better)"),
    Height_comparison %>%
      mutate(Fit = "Cover x height", Panel = recode(Term, Cover = "Canopy cover (WVCC)", Height = "Woody canopy height (WVCH)", Interaction = "Cover x height interaction"))
  ) %>%
    mutate(Panel = factor(Panel, c("WAIC minus baseline WAIC (lower is better)", "Canopy cover (WVCC)", "Woody canopy height (WVCH)", "Cover x height interaction"))) %>%
    filter(!is.na(Median)) %>%
    ggplot(aes(Radius_m, Median, colour = Fit)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_linerange(aes(ymin = Lo, ymax = Hi), na.rm = TRUE) +
    geom_line() +
    geom_point() +
    facet_wrap(~Panel, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = Height_fit_radii_m, labels = radius_label) +
    scale_colour_manual(values = c(`Cover x height` = "black", `Cover only` = "grey60")) +
    labs(x = paste0("Ring outer radius (inner ", Ring_inner_m, " m excluded)"), y = "Community-level mean (logit, per 1 SD), 95% CrI", colour = NULL) +
    theme(legend.position = "top")
  ggsave(file.path(Fig_dir, "Scale_of_effect_height_radii.png"), height_radii_plot, width = 7, height = 9, bg = "white")
  print(height_radii_plot)

  ### Interaction: community-level mean partial effect of ring cover on the logit scale (relative to the ring's mean cover at its mean height), over the central 90% of observed cover, at the 10th, 50th and 90th percentiles of woody canopy height; covariate scaling as in the fits (over the surveyed rows of each radius, rings without canopy left out of the height scaling)
  Height_scaling <- Rows_height %>%
    filter(Radius_m %in% as.numeric(names(Height_summaries))) %>%
    summarize(.by = Radius_m, Cover_mean = mean(Canopy_cover_check), Cover_sd = sd(Canopy_cover_check), Cover_lo = quantile(Canopy_cover_check, 0.05), Cover_hi = quantile(Canopy_cover_check, 0.95),
              Height_mean = mean(Canopy_height_woody_ring, na.rm = TRUE), Height_sd = sd(Canopy_height_woody_ring, na.rm = TRUE),
              Height_q = list(quantile(Canopy_height_woody_ring, c(0.1, 0.5, 0.9), na.rm = TRUE)))
  Interaction_curves <- Height_scaling %>%
    pmap(\(Radius_m, Cover_mean, Cover_sd, Cover_lo, Cover_hi, Height_mean, Height_sd, Height_q) {
      b <- Height_summaries[[as.character(Radius_m)]]$draws$beta.comm
      b_cover <- b[, term_name("Cover_h", Radius_m)]; b_height <- b[, term_name("Height_h", Radius_m)]; b_int <- b[, term_name("Interaction", Radius_m)]
      expand_grid(Cover = seq(Cover_lo, Cover_hi, length.out = 50), Height_pctl = names(Height_q)) %>%
        mutate(Height_m = Height_q[Height_pctl], draws = map2(Cover, Height_m, \(cv, ht) {
          cz <- (cv - Cover_mean) / Cover_sd; hz <- (ht - Height_mean) / Height_sd
          b_cover * cz + b_height * hz + b_int * cz * hz
        }), Radius_m = Radius_m)
    }) %>%
    list_rbind() %>%
    mutate(Median = map_dbl(draws, median), Lo = map_dbl(draws, quantile, 0.025), Hi = map_dbl(draws, quantile, 0.975),
           Height_label = paste0(sub("%", "th", Height_pctl), " percentile"),
           Panel = factor(radius_label(Radius_m), radius_label(sort(unique(Radius_m))))) %>%
    select(-draws)
  ## Cover is a 0-1 share in the height table; shown as percent
  height_interaction_plot <- ggplot(Interaction_curves, aes(Cover * 100, Median, colour = Height_label, fill = Height_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_text(data = distinct(Interaction_curves, Panel, Height_label, Height_m) %>% mutate(lab = sprintf("%.1f m", Height_m)) %>% summarize(.by = Panel, lab = paste(lab, collapse = " / ")), aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 3) +
    facet_wrap(~Panel, scales = "free") +
    scale_colour_viridis_d(end = 0.85, aesthetics = c("colour", "fill")) +
    labs(x = "Canopy cover in the ring (%)", y = "Community-level mean partial effect\n(logit, relative to mean cover and height)", colour = "Woody canopy height", fill = "Woody canopy height") +
    theme(legend.position = "top")
  ggsave(file.path(Fig_dir, "Scale_of_effect_height_interaction.png"), height_interaction_plot, width = 8, height = 7, bg = "white")
  print(height_interaction_plot)

  print(Height_comparison, n = Inf, width = Inf)
}
