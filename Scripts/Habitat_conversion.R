## PhD birds in silvopastoral landscapes ##
# Habitat-conversion scenarios (pasture to silvopasture and back) at the surveyed locations, predicted from the habitat-conversion model (habitat by 1.5 km canopy cover, by silvopasture within 300 m and by rainfall, with the cover-by-height adjustment and forest cover within 300 m, so the habitat-specific ssp slopes hold forest fixed as the scenarios do).
# Each scenario changes a location's habitat, its silvopasture cover within 300 m (ssp) and, in the buffer scenarios, its canopy cover within 1.5 km, together with every interaction column built from them; all other covariates (woody canopy height within 1.5 km included) keep their observed values, ssp is kept between 0 and the highest value observed at the point's habitat after conversion (each habitat has its own ssp slope, estimated only over that habitat's range), and canopy cover within its observed range.
# With the canopy-weighted silvopasture fit (the default), pasture converted to silvopasture adds canopy-weighted silvopasture at the canopy density (canopy-weighted / area-based silvopasture) of a randomly drawn farm, drawn independently for each location and posterior draw, so the intervals carry both model uncertainty and the observed variation among farms in the kind of silvopasture they plant; farms, not 300 m buffers, are the sampling unit, as the farmer decides the silvopasture type and overlapping buffers would count shared stands twice. Silvopasture converted to pasture removes the canopy-weighted silvopasture actually mapped there.
# Counterfactual occupancy comes from the fitted occupancy probability of each species at each surveyed location-period and draw (which already carries the location's spatial effect) by adding the change in the linear predictor on the logit scale; the latent occupancy states are not used, as they are fixed at 1 wherever a species was detected and so cannot respond to a scenario.
# Change in expected richness = sum over the modelled species of counterfactual minus fitted occupancy, per draw; summaries are the mean over the changed location-periods per draw, then the posterior median and CrI.
# Overall change in landscape richness, ignoring the sampling design: expected number of species across one average location-period of each habitat, sum over species of 1 - product over habitats of (1 - the habitat's mean occupancy), so species shared between habitats count once and each habitat weighs the same; assumes a species' presence is independent among the three habitats.
# For the locations whose habitat switches, the share of the gap to the real locations of the new habitat that conversion closes, per draw: in expected richness, (after - before) / (target - before); in composition, 1 - Bray-Curtis(after, target) / Bray-Curtis(before, target), where each community is the mean occupancy of every species over the group's location-periods (Bray-Curtis on occupancy probabilities, so it responds to richness as well as turnover).
# Composition is drawn as an NMDS of the real location-periods, held fixed, with each scenario's converted location-periods added to it (vegan::MDSaddpoints), so the simulated communities never shape the ordination of the real ones.

# Libraries ---------------------------------------------------------------
library(tidyverse)
ggplot2::theme_set(cowplot::theme_cowplot())

# Configuration -----------------------------------------------------------
## Fit driving the scenarios, from the command line: 'f300_cw' (default; the habitat-conversion model with forest cover within 300 m and canopy-weighted silvopasture within 300 m), 'f300' (the same with area-based silvopasture) or 'f300_rip' (area-based, plus a riparian-forest indicator, which keeps its observed value in every scenario, as forest points are never converted); a variant other than the default adds its tag to every output CSV and writes its figures to its own folder
Fit_variant <- c(commandArgs(trailingOnly = TRUE), "f300_cw")[1]
Occ_mod_path <- c(
  f300_cw = "Derived/Rdata/Model_output/Temporal/Scale_of_effect/Occ_mod_98spp_soe_full1500_f300_cw_2026-11-23.rds",
  f300 = "Derived/Rdata/Model_output/Temporal/Scale_of_effect/Occ_mod_98spp_soe_full1500_f300_2026-11-22.rds",
  f300_rip = "Derived/Rdata/Model_output/Temporal/Scale_of_effect/Occ_mod_98spp_soe_full1500_f300_rip_2026-11-22.rds"
)[[Fit_variant]]
Out_suffix <- if (Fit_variant == "f300_cw") "" else paste0("_", Fit_variant)
# The fit's silvopasture-within-300 m column: canopy-weighted (Ssp_cw) or area-based (ssp)
Ssp_col <- c(f300_cw = "Ssp_cw", f300 = "ssp", f300_rip = "ssp")[[Fit_variant]]
Ssp_label <- c(Ssp_cw = "Silvopasture canopy within 300 m (%)", ssp = "Silvopasture within 300 m (%)")[[Ssp_col]]
# Canopy-weighted silvopasture within 300 m and within the 50 m disc of each survey, and each point count's farm (Id_group without its institution prefix), written by the silvopasture-canopy script
Ssp_canopy_csv <- "Derived/Excels/Ssp_canopy_by_survey.csv"
Farm_csv <- "Derived/Excels/Isolated_trees_by_location.csv"
# Farms whose mean area-based silvopasture within 300 m is below this (percent) give no canopy density, as their ratio rests on a sliver of silvopasture
Min_farm_ssp_pct <- 2
Density_seed <- 20261124
# Stacked list the multi-season list was built from (one row per surveyed location-period, with the survey ID that links to the landscape metrics) and its scaling table
msOcc_l_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm.rds"
# Multi-season list the fit was built from, for the location and period order of the fit's arrays (which carry no names)
tmsOcc_l_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_98_hab_fix.rds"
# Land cover within 300 m of each point (percent forest, intensive pasture, silvopasture, other) per visit
Lsm_path <- "../Farm_diversity/Derived/Excels/Event_covs_lsm.csv"
# The same land cover at radii from 25 to 300 m, one workbook per digitization (the 300 m values link each surveyed row to its workbook row, as the survey table carries the 300 m values only)
Lsm_radii_dir <- "../Farm_diversity/Derived/Excels/Lsm"
Lsm_radii_pattern <- "^Lsm_df_.*_10\\.29\\.26\\.xlsx$"
# Canopy cover within 1.5 km (ring from 50 m out, percent) of each surveyed row, written by the scale-of-effect script, and the radius the model uses
Fit_covs_path <- "Derived/Excels/Scale_of_effect_fit_covs.csv"
Cover_radius_m <- 1500
# Point coordinates and farm, and the canopy rasters (farm x year clips, else the national map) used for the ring covers
Site_covs_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Site_covs.csv"
Clip_dir <- "../Farm_diversity/Data/Geospatial/canopy_rasters"
National_dir <- "../Geospatial_data/Environmental"
Raster_years <- c(2013, 2014, 2016, 2017, 2019, 2022, 2024)
Utm_crs <- "EPSG:32618"

# Converting the pasture (or silvopasture) inside a point's own Disc_radius_m disc changes ssp, the percent of the Buffer_radius_m disc, by that class's percent of the small disc times the small disc's share of the buffer
Disc_radius_m <- 50
Buffer_radius_m <- 300
Disc_share_pct <- 100 * (Disc_radius_m / Buffer_radius_m)^2

## Scenarios: 'disc' converts each point's own 50 m disc (pasture points become silvopasture or the reverse); 'buffer' converts all pasture (or all silvopasture) within every point's 300 m buffer, so forest points' ssp changes too
# Converting land changes canopy cover by the difference in mean canopy cover between the 50 m discs of silvopasture and pasture points whose disc is at least Pure_disc_pct of their own class; the 1.5 km ring starts at 50 m, so a point's own disc never enters its ring and the disc scenarios leave canopy cover unchanged
# The buffer scenarios always carry the resulting change in canopy cover within 1.5 km (user, 2026-11-24: the held-cover versions changed the results little)
Scenarios <- tribble(
  ~Scenario,                       ~Direction,   ~Extent,  ~Change_cover,
  "disc_to_ssp",                   "to_ssp",     "disc",   FALSE,
  "disc_to_pasture",               "to_pasture", "disc",   FALSE,
  "buffer_to_ssp",                 "to_ssp",     "buffer", TRUE,
  "buffer_to_pasture",             "to_pasture", "buffer", TRUE
)
Pure_disc_pct <- 90
ci_level <- 0.95
## Composition ordination: posterior draws averaged into the Bray-Curtis dissimilarities, NMDS dimensions and random starts (as in the NMDS script)
Ord_n_draws <- 50
Ord_k <- 3
Ord_try <- 20
Ord_seed <- 123

Disc_canopy_csv <- "Derived/Excels/Canopy_disc50_by_location.csv"
Out_csv <- paste0("Derived/Excels/Habitat_conversion_summary", Out_suffix, ".csv")
Out_species_csv <- paste0("Derived/Excels/Habitat_conversion_species", Out_suffix, ".csv")
Out_locations_csv <- paste0("Derived/Excels/Habitat_conversion_locations", Out_suffix, ".csv")
Out_landscape_csv <- paste0("Derived/Excels/Habitat_conversion_landscape", Out_suffix, ".csv")
Fig_dir <- c(f300_cw = "Figures/Habitat_conversion", f300 = "Figures/Habitat_conversion/Area_ssp", f300_rip = "Figures/Riparian/Habitat_conversion")[[Fit_variant]]

level_labels <- c(Pastizales = "Pasture", Ssp = "Silvopasture", Bosque = "Forest")
level_colors <- c(Pasture = "#eda100", Silvopasture = "#6250d6", Forest = "#008300")
Scenario_labels <- c(
  disc_to_ssp = "Pasture to silvopasture,\npoint count circle (50 m)", disc_to_pasture = "Silvopasture to pasture,\npoint count circle (50 m)",
  buffer_to_ssp = "All pasture to silvopasture,\nwithin 300 m", buffer_to_pasture = "All silvopasture to pasture,\nwithin 300 m"
)

# Functions ---------------------------------------------------------------
# New habitat and canopy cover (percent) of each surveyed location-period under one scenario, and its new silvopasture within 300 m as Ssp_base + density * Ssp_gain, where Ssp_gain is the area of pasture converted to silvopasture (percent of the 300 m buffer) and density the canopy-weighted silvopasture each unit of it adds (1 for area-based silvopasture); cells has one row per location-period with Habitat, the fitted silvopasture Ssp_fit (percent of the 300 m buffer) and Ssp_fit_50 (percent of the 50 m disc), intpast and intpast_50 (pasture, percent of the buffer and disc), Cover_pct and the canopy-cover change of each direction (dCover_to_ssp, dCover_to_pasture); ssp_max has the highest observed Ssp_fit of each habitat (Habitat, ssp_max), cover_range the observed range of Cover_pct
scenario_covs <- function(cells, direction, extent, change_cover, disc_share_pct, ssp_max, cover_range) {
  from <- c(to_ssp = "Pastizales", to_pasture = "Ssp")[[direction]]
  to <- c(to_ssp = "Ssp", to_pasture = "Pastizales")[[direction]]
  convert <- cells$Habitat == from
  ssp_base <- if (extent == "disc" && direction == "to_pasture") {
    cells$Ssp_fit - convert * disc_share_pct * cells$Ssp_fit_50 / 100
  } else if (direction == "to_pasture") {
    rep(0, nrow(cells))
  } else {
    cells$Ssp_fit
  }
  ssp_gain <- if (direction == "to_pasture") rep(0, nrow(cells)) else if (extent == "disc") convert * disc_share_pct * cells$intpast_50 / 100 else cells$intpast
  cover_new <- if (change_cover) cells$Cover_pct + cells[[paste0("dCover_", direction)]] else cells$Cover_pct
  mutate(cells, Habitat_new = if_else(convert, to, Habitat), Ssp_base = ssp_base, Ssp_gain = ssp_gain, cover_raw = cover_new) %>%
    left_join(ssp_max, by = c(Habitat_new = "Habitat"), relationship = "many-to-one") %>%
    mutate(Cover_new = pmin(pmax(cover_raw, cover_range[1]), cover_range[2]))
}

# Bray-Curtis dissimilarity between matching rows of two draws x species matrices of occupancy probabilities
bray_curtis_rows <- function(comm_a, comm_b) rowSums(abs(comm_a - comm_b)) / rowSums(comm_a + comm_b)

# Load data ---------------------------------------------------------------
Occ_mod <- readRDS(Occ_mod_path)
X <- Occ_mod$X
psi <- Occ_mod$psi.samples
beta <- Occ_mod$beta.samples
sp_names <- Occ_mod$sp.names
rm(Occ_mod)
invisible(gc())
Site_levels <- dimnames(readRDS(tmsOcc_l_path)$y)[[2]]
Time_levels <- dimnames(readRDS(tmsOcc_l_path)$y)[[3]]
stopifnot(length(Site_levels) == dim(X)[1], length(Time_levels) == dim(X)[2])

## Surveyed location-periods of the multi-season fit, linked to their survey ID through the stacked list (site name = location code, period = year group and season); the survey year is the latest within the location-period, snapped to the nearest canopy raster year, as for the ring covers
Scaling <- read_csv(str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv"), show_col_types = FALSE) %>%
  filter(Process == "occupancy")
Lsm_rows <- read_csv(Lsm_path, col_types = cols(Id_survey = "c", Year_grp = "c", Season = "c")) %>%
  slice_max(by = c(Id_survey, Year_grp, Season), order_by = Year, with_ties = FALSE) %>%
  select(Id_survey, Id_survey_no_dc, Year_grp, Season, Year, forest, ssp, intpast)
Cells <- readRDS(msOcc_l_path)$occ.covs %>%
  select(Id_survey, Year_grp, Season, Site_code = Id_survey_no_dc, Habitat) %>%
  mutate(across(everything(), as.character), Site = match(Site_code, Site_levels), Period = match(paste(Year_grp, Season, sep = "."), Time_levels)) %>%
  left_join(Lsm_rows, by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
  left_join(read_csv(Fit_covs_path, col_types = cols(Id_survey = "c", Year_grp = "c", Season = "c")) %>% filter(Radius_m == Cover_radius_m) %>% select(Id_survey, Year_grp, Season, Cover_pct), by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
  mutate(Raster_year = map_dbl(Year, \(year) Raster_years[which.min(abs(Raster_years - year))]))
stopifnot(!anyNA(Cells$Site), !anyNA(Cells$Period), !anyNA(Cells$intpast), !anyNA(Cells$Cover_pct), !anyDuplicated(Cells[c("Site", "Period")]))
## Land cover in each row's 50 m disc, from the workbook row whose 300 m values match the survey table's (the closest, where two digitizations of a location agree to within rounding)
Lsm_radii <- map(list.files(Lsm_radii_dir, Lsm_radii_pattern, full.names = TRUE), readxl::read_excel) %>%
  list_rbind() %>%
  rename(Src_id = Id_survey) %>%
  left_join(distinct(Lsm_rows, Src_id = Id_survey, Id_survey_no_dc), by = "Src_id", relationship = "many-to-one")
Disc_cover <- Cells %>%
  select(Id_survey, Year_grp, Season, Id_survey_no_dc, forest, intpast, ssp) %>%
  inner_join(filter(Lsm_radii, buffer == Buffer_radius_m) %>% select(Id_survey_no_dc, Src_id, lc_file, data_year, f3 = forest, p3 = intpast, s3 = ssp), by = "Id_survey_no_dc", relationship = "many-to-many") %>%
  mutate(Diff = abs(f3 - forest) + abs(p3 - intpast) + abs(s3 - ssp)) %>%
  filter(Diff < 0.15) %>%
  slice_min(Diff, by = c(Id_survey, Year_grp, Season), with_ties = FALSE) %>%
  left_join(filter(Lsm_radii, buffer == Disc_radius_m) %>% select(Src_id, lc_file, data_year, forest_50 = forest, intpast_50 = intpast, ssp_50 = ssp), by = c("Src_id", "lc_file", "data_year"), relationship = "many-to-one") %>%
  select(Id_survey, Year_grp, Season, forest_50, intpast_50, ssp_50)
## Canopy-weighted silvopasture within 300 m and within the 50 m disc, and the farm of each row; Ssp_fit and Ssp_fit_50 are the silvopasture measure the fit uses
Farms <- read_csv(Farm_csv, show_col_types = FALSE) %>%
  mutate(Farm = str_remove(Id_group, "^[A-Z]+-")) %>%
  select(Id_survey, Farm)
Ssp_canopy <- read_csv(Ssp_canopy_csv, col_types = cols(Id_survey = "c", Year_grp = "c"))
Cells <- left_join(Cells, Disc_cover, by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
  left_join(select(Ssp_canopy, Id_survey, Year_grp, Ssp_cw, Ssp_cw_50), by = c("Id_survey", "Year_grp"), relationship = "many-to-one") %>%
  left_join(Farms, by = "Id_survey", relationship = "many-to-one") %>%
  mutate(Ssp_fit = .data[[Ssp_col]], Ssp_fit_50 = if (Ssp_col == "Ssp_cw") Ssp_cw_50 else ssp_50)
stopifnot(!anyNA(Cells$intpast_50), !anyNA(Cells$ssp_50), !anyNA(Cells$Ssp_cw_50), !anyNA(Cells$Farm))
Cover_scaling <- tibble(Mean = mean(Cells$Cover_pct), SD = sd(Cells$Cover_pct))
## Area-based ssp was scaled with the stacked list, canopy-weighted ssp over the surveyed location-periods when the fit was built
Ssp_scaling <- if (Ssp_col == "ssp") filter(Scaling, Covariate == "ssp") else tibble(Mean = mean(Cells$Ssp_fit), SD = sd(Cells$Ssp_fit))
X_obs <- t(apply(X, 3, \(m) m[cbind(Cells$Site, Cells$Period)]))
## The landscape metrics and ring covers must reproduce the fitted ssp, canopy cover and habitat, so the scenarios start from the values the model saw (ring cover is stored to 3 decimals)
stopifnot(max(abs(X_obs[Ssp_col, ] * Ssp_scaling$SD + Ssp_scaling$Mean - Cells$Ssp_fit)) < 1e-6, max(abs(X_obs["Cover_h1500", ] * Cover_scaling$SD + Cover_scaling$Mean - Cells$Cover_pct)) < 1e-3)
stopifnot(all(X_obs["HabitatSsp", ] == (Cells$Habitat == "Ssp")), all(X_obs["HabitatBosque", ] == (Cells$Habitat == "Bosque")))
X_obs <- t(X_obs)

## Fitted occupancy on the logit scale, draws x species x surveyed location-periods (bounded away from 0 and 1 so the logit is finite)
cell_index <- Cells$Site + (Cells$Period - 1) * dim(X)[1]
Logit_obs <- qlogis(pmin(pmax(array(psi, c(dim(psi)[1:2], prod(dim(psi)[3:4])))[, , cell_index], 1e-12), 1 - 1e-12))
Psi_obs <- plogis(Logit_obs)
rm(psi)
invisible(gc())

# Canopy cover change from conversion ----
## Mean canopy cover in the 50 m disc of every pasture and silvopasture location in its survey years, from the same rasters as the ring covers; run once, then read from the CSV
Locs <- read_csv(Site_covs_path, show_col_types = FALSE) %>%
  distinct(Id_survey_no_dc, Id_scr, Long, Lat)
if (!file.exists(Disc_canopy_csv)) {
  Cells %>%
    filter(Habitat %in% c("Pastizales", "Ssp")) %>%
    distinct(Id_survey_no_dc, Raster_year) %>%
    left_join(Locs, by = "Id_survey_no_dc", relationship = "many-to-one") %>%
    group_split(Id_scr, Raster_year) %>%
    map(\(grp) {
      pts <- terra::vect(grp, geom = c("Long", "Lat"), crs = "EPSG:4326")
      clip_path <- file.path(Clip_dir, sprintf("farm_%s_%d.tif", grp$Id_scr[1], grp$Raster_year[1]))
      canopy_r <- if (file.exists(clip_path)) terra::rast(clip_path) else terra::crop(terra::rast(file.path(National_dir, sprintf("WVCC/Colombia_WVCC_%d_50km.tif", grp$Raster_year[1]))), terra::buffer(pts, 500), snap = "out")
      mutate(select(grp, Id_survey_no_dc, Raster_year), Canopy_disc = exactextractr::exact_extract(canopy_r, sf::st_as_sf(terra::buffer(pts, width = Disc_radius_m)), "mean", progress = FALSE))
    }) %>%
    list_rbind() %>%
    mutate(Canopy_disc = round(Canopy_disc, 3)) %>%
    arrange(Id_survey_no_dc, Raster_year) %>%
    write_csv(Disc_canopy_csv)
}
Disc_canopy <- read_csv(Disc_canopy_csv, show_col_types = FALSE)
stopifnot(!anyNA(Disc_canopy$Canopy_disc), max(Disc_canopy$Canopy_disc) <= 100)
## Canopy cover gained where pasture becomes silvopasture (percentage points; the rasters store percent), over the surveyed location-periods of each habitat
Disc_by_habitat <- Cells %>%
  filter((Habitat == "Ssp" & ssp_50 >= Pure_disc_pct) | (Habitat == "Pastizales" & intpast_50 >= Pure_disc_pct)) %>%
  inner_join(Disc_canopy, by = c("Id_survey_no_dc", "Raster_year"), relationship = "many-to-one") %>%
  summarise(Canopy_disc = mean(Canopy_disc), N = n(), .by = Habitat)
Canopy_gain <- with(Disc_by_habitat, Canopy_disc[Habitat == "Ssp"] - Canopy_disc[Habitat == "Pastizales"])

## Area of the 1.5 km ring (50 m out) covered by the 300 m buffers of all surveyed locations, and the shares of pasture and silvopasture in it (each buffer's share averaged over its surveyed periods, buffers weighted by their overlap with the ring)
Loc_shares <- Cells %>%
  summarise(intpast = mean(intpast), ssp = mean(ssp), .by = c(Site, Id_survey_no_dc)) %>%
  left_join(Locs, by = "Id_survey_no_dc", relationship = "many-to-one") %>%
  arrange(Site) %>%
  sf::st_as_sf(coords = c("Long", "Lat"), crs = "EPSG:4326") %>%
  sf::st_transform(Utm_crs)
Buffers <- sf::st_geometry(sf::st_buffer(Loc_shares, Buffer_radius_m))
Rings <- sf::st_sfc(map2(sf::st_geometry(sf::st_buffer(Loc_shares, Cover_radius_m)), sf::st_geometry(sf::st_buffer(Loc_shares, Disc_radius_m)), sf::st_difference), crs = Utm_crs)
Ring_overlap <- imap(sf::st_intersects(Rings, Buffers), \(nb, i) {
  overlap <- map_dbl(nb, \(j) sum(as.numeric(sf::st_area(sf::st_intersection(Buffers[j], Rings[i])))))
  tibble(Site = Loc_shares$Site[i], Ring_area = as.numeric(sf::st_area(Rings[i])), Union_area = as.numeric(sf::st_area(sf::st_intersection(sf::st_union(Buffers[nb]), Rings[i]))),
    Pasture_share = weighted.mean(Loc_shares$intpast[nb], overlap) / 100, Ssp_share = weighted.mean(Loc_shares$ssp[nb], overlap) / 100)
}) %>% list_rbind()
## Change in ring canopy cover (percentage points) when every buffer's pasture becomes silvopasture, or every buffer's silvopasture becomes pasture
Cells <- Cells %>%
  left_join(Ring_overlap, by = "Site", relationship = "many-to-one") %>%
  mutate(dCover_to_ssp = Union_area * Pasture_share * Canopy_gain / Ring_area, dCover_to_pasture = -Union_area * Ssp_share * Canopy_gain / Ring_area)

# Scenarios ---------------------------------------------------------------
## Highest observed silvopasture within 300 m of each habitat, in the fit's measure
Ssp_max <- summarise(Cells, ssp_max = max(Ssp_fit), .by = Habitat)
## Canopy density of each farm's silvopasture (canopy-weighted / area-based silvopasture, summed over its point counts' buffers and year-groups), and the density each location receives in each posterior draw, one farm drawn per location and draw (1 throughout for area-based silvopasture)
Farm_density <- Ssp_canopy %>%
  inner_join(Farms, by = "Id_survey", relationship = "many-to-one") %>%
  summarise(Density = sum(Ssp_cw) / sum(Ssp_area_poly), Ssp_pct = mean(Ssp_area_poly), .by = Farm) %>%
  filter(Ssp_pct >= Min_farm_ssp_pct) %>%
  arrange(Farm)
set.seed(Density_seed)
Density <- if (Ssp_col == "Ssp_cw") {
  matrix(sample(Farm_density$Density, nrow(beta) * length(Site_levels), replace = TRUE), nrow(beta))[, Cells$Site]
} else {
  matrix(1, nrow(beta), nrow(Cells))
}
probs <- c((1 - ci_level) / 2, 0.5, 1 - (1 - ci_level) / 2)
Prec_scaling <- filter(Scaling, Covariate == "Tot_prec")
## Expected richness and mean community (mean occupancy of each species) over the real location-periods of each habitat, per draw: the targets converted locations are compared with
Habitat_cells <- split(seq_len(nrow(Cells)), unname(level_labels[Cells$Habitat]))
Rich_by_habitat <- map(Habitat_cells, \(idx) rowMeans(apply(Psi_obs[, , idx, drop = FALSE], c(1, 3), sum)))
Comm_by_habitat <- map(Habitat_cells, \(idx) rowMeans(Psi_obs[, , idx, drop = FALSE], dims = 2))
Ord_draws <- round(seq(1, dim(Psi_obs)[1], length.out = Ord_n_draws))
Scenario_results <- pmap(Scenarios, \(Scenario, Direction, Extent, Change_cover) {
  New <- scenario_covs(Cells, Direction, Extent, Change_cover, Disc_share_pct, Ssp_max, range(Cells$Cover_pct))
  From_label <- c(to_ssp = "Pasture", to_pasture = "Silvopasture")[[Direction]]
  To_label <- c(to_ssp = "Silvopasture", to_pasture = "Pasture")[[Direction]]
  ## New silvopasture within 300 m per draw and location-period, kept between 0 and the highest value observed at the new habitat; a location-period changes if its habitat, canopy cover or silvopasture changes at a density of 1 (the highest a farm can have)
  ssp_raw <- sweep(Density * rep(New$Ssp_gain, each = nrow(Density)), 2, New$Ssp_base, `+`)
  ssp_new <- pmin(pmax(ssp_raw, 0), rep(New$ssp_max, each = nrow(Density)))
  ssp_new_1 <- pmin(pmax(New$Ssp_base + New$Ssp_gain, 0), New$ssp_max)
  ## Counterfactual design matrix: habitat dummies and scaled canopy cover replaced, every interaction column recomputed as the product of its parts; silvopasture varies by draw, so its column is set to 1 here and the columns containing it (Ssp_cols) are multiplied by each draw's scaled silvopasture below
  X_new <- X_obs
  X_new[, "HabitatSsp"] <- as.numeric(New$Habitat_new == "Ssp")
  X_new[, "HabitatBosque"] <- as.numeric(New$Habitat_new == "Bosque")
  X_new[, Ssp_col] <- 1
  X_new[, "Cover_h1500"] <- (New$Cover_new - Cover_scaling$Mean) / Cover_scaling$SD
  for (col in str_subset(colnames(X_new), ":")) X_new[, col] <- apply(X_new[, strsplit(col, ":", fixed = TRUE)[[1]], drop = FALSE], 1, prod)
  Ssp_cols <- colnames(X_new)[map_lgl(strsplit(colnames(X_new), ":", fixed = TRUE), \(parts) Ssp_col %in% parts)]
  Other_cols <- setdiff(colnames(X_new), Ssp_cols)
  changed <- which(rowSums(abs(X_new[, Other_cols, drop = FALSE] - X_obs[, Other_cols, drop = FALSE])) > 1e-12 | abs(ssp_new_1 - Cells$Ssp_fit) > 1e-12)
  dX_other <- X_new[changed, Other_cols, drop = FALSE] - X_obs[changed, Other_cols, drop = FALSE]
  ssp_scaled <- (ssp_new[, changed, drop = FALSE] - Ssp_scaling$Mean) / Ssp_scaling$SD

  ## Change in occupancy per draw, changed location-period and species; species coefficients are named "<coefficient>-<species>"
  dPsi <- vapply(seq_along(sp_names), \(s) {
    b <- beta[, paste(colnames(X_obs), sp_names[s], sep = "-"), drop = FALSE]
    colnames(b) <- colnames(X_obs)
    d_eta <- b[, Other_cols, drop = FALSE] %*% t(dX_other) + ssp_scaled * (b[, Ssp_cols, drop = FALSE] %*% t(X_new[changed, Ssp_cols, drop = FALSE])) - b[, Ssp_cols, drop = FALSE] %*% t(X_obs[changed, Ssp_cols, drop = FALSE])
    plogis(Logit_obs[, s, changed, drop = FALSE][, 1, ] + d_eta) - Psi_obs[, s, changed, drop = FALSE][, 1, ]
  }, matrix(0, nrow(beta), length(changed)))
  dRich <- apply(dPsi, c(1, 2), sum)
  Rich_obs <- apply(Psi_obs[, , changed, drop = FALSE], c(1, 3), sum)

  ## Mean over the changed location-periods per draw, for all of them and by their original habitat
  habitat_from <- unname(level_labels[Cells$Habitat[changed]])
  groups <- c(list(All = seq_along(changed)), split(seq_along(changed), habitat_from))
  summary_df <- imap(groups, \(idx, group) {
    d <- rowMeans(dRich[, idx, drop = FALSE])
    q <- quantile(d, probs)
    rich_before <- rowMeans(Rich_obs[, idx, drop = FALSE])
    q_after <- quantile(rich_before + d, probs)
    ## Share of the gap to the new habitat closed, only for the group whose habitat switched
    if (group == From_label) {
      q_gap <- quantile(d / (Rich_by_habitat[[To_label]] - rich_before), probs)
      comm_before <- rowMeans(Psi_obs[, , changed[idx], drop = FALSE], dims = 2)
      comm_after <- comm_before + apply(dPsi[, idx, , drop = FALSE], c(1, 3), mean)
      bc_before <- bray_curtis_rows(comm_before, Comm_by_habitat[[To_label]])
      bc_after <- bray_curtis_rows(comm_after, Comm_by_habitat[[To_label]])
      q_comp <- quantile(1 - bc_after / bc_before, probs)
      target <- c(median(Rich_by_habitat[[To_label]]), median(bc_before), median(bc_after))
    } else {
      q_gap <- q_comp <- rep(NA_real_, 3)
      target <- rep(NA_real_, 3)
    }
    tibble(Scenario = Scenario, Habitat_from = group, N_cells = length(idx), N_locations = n_distinct(Cells$Site[changed][idx]),
      Rich_obs = median(rich_before), dRich_low = q[[1]], dRich_median = q[[2]], dRich_high = q[[3]], P_gt0 = mean(d > 0),
      dSsp_mean = mean(sweep(ssp_new[, changed[idx], drop = FALSE], 2, Cells$Ssp_fit[changed][idx])), dCover_mean = mean(New$Cover_new[changed][idx] - Cells$Cover_pct[changed][idx]), N_ssp_capped = mean(rowSums(ssp_new[, changed[idx], drop = FALSE] != ssp_raw[, changed[idx], drop = FALSE])),
      Rich_after_low = q_after[[1]], Rich_after_median = q_after[[2]], Rich_after_high = q_after[[3]], Rich_target = target[1],
      Gap_closed_low = q_gap[[1]], Gap_closed_median = q_gap[[2]], Gap_closed_high = q_gap[[3]],
      BC_before = target[2], BC_after = target[3], Comp_closed_low = q_comp[[1]], Comp_closed_median = q_comp[[2]], Comp_closed_high = q_comp[[3]])
  }) %>% list_rbind()
  ## Species: mean change in occupancy over the changed location-periods of each original habitat, per draw
  species_df <- imap(groups[-1], \(idx, group) {
    sp_draws <- apply(dPsi[, idx, , drop = FALSE], c(1, 3), mean)
    tibble(Scenario = Scenario, Habitat_from = group, Species = sp_names, dPsi_low = apply(sp_draws, 2, quantile, probs[1]), dPsi_median = apply(sp_draws, 2, median), dPsi_high = apply(sp_draws, 2, quantile, probs[3]), P_gt0 = colMeans(sp_draws > 0))
  }) %>% list_rbind()
  ## Locations: change averaged over each changed location's surveyed periods per draw, with the location's observed covariates (original units, averaged over its periods)
  site_draws <- t(rowsum(t(dRich), Cells$Site[changed])) / rep(as.vector(table(Cells$Site[changed])), each = nrow(dRich))
  location_df <- tibble(Site = as.integer(colnames(site_draws)), dRich_low = apply(site_draws, 2, quantile, probs[1]), dRich_median = apply(site_draws, 2, median), dRich_high = apply(site_draws, 2, quantile, probs[3])) %>%
    left_join(Cells[changed, ] %>%
      mutate(Tot_prec = X_obs[changed, "Tot_prec"] * Prec_scaling$SD + Prec_scaling$Mean) %>%
      summarise(Habitat_from = unname(level_labels[Habitat[1]]), Cover_pct = mean(Cover_pct), !!Ssp_col := mean(Ssp_fit), Tot_prec = mean(Tot_prec), .by = Site), by = "Site") %>%
    mutate(Scenario = Scenario, .before = 1)
  ## Landscape richness ignoring the sampling design: expected number of species across one average location-period of each original habitat (the union of the three habitats' mean communities, so a species shared between habitats counts once, and each habitat weighs the same whatever its number of surveyed location-periods), before and after, per draw
  comm_after <- imap(Comm_by_habitat, \(comm, habitat) {
    idx <- which(habitat_from == habitat)
    if (!length(idx)) return(comm)
    comm + apply(dPsi[, idx, , drop = FALSE], c(1, 3), sum) / length(Habitat_cells[[habitat]])
  })
  union_richness <- \(comms) rowSums(1 - reduce(map(comms, \(comm) 1 - comm), `*`))
  rich_landscape <- union_richness(Comm_by_habitat)
  d_landscape <- union_richness(comm_after) - rich_landscape
  q_landscape <- quantile(d_landscape, probs)
  landscape_df <- tibble(Scenario = Scenario, Rich_landscape_before = median(rich_landscape), dRich_landscape_low = q_landscape[[1]], dRich_landscape_median = q_landscape[[2]], dRich_landscape_high = q_landscape[[3]], P_gt0 = mean(d_landscape > 0))
  ## Counterfactual occupancy of the changed location-periods at the ordination draws, draws x species x changed location-periods
  psi_new <- Psi_obs[Ord_draws, , changed, drop = FALSE] + aperm(dPsi[Ord_draws, , , drop = FALSE], c(1, 3, 2))
  list(summary = summary_df, species = species_df, locations = location_df, landscape = landscape_df, changed = changed, psi_new = psi_new)
})
names(Scenario_results) <- Scenarios$Scenario

# Outputs -----------------------------------------------------------------
Summary <- map(Scenario_results, "summary") %>%
  list_rbind() %>%
  mutate(across(c(Rich_obs, starts_with("dRich"), dSsp_mean, dCover_mean, starts_with("Rich_after"), Rich_target, starts_with("Gap_closed"), starts_with("BC_"), starts_with("Comp_closed")), \(x) round(x, 3)), P_gt0 = round(P_gt0, 3))
write_csv(Summary, Out_csv)
Species <- map(Scenario_results, "species") %>%
  list_rbind() %>%
  mutate(across(starts_with("dPsi"), \(x) round(x, 4)), P_gt0 = round(P_gt0, 3)) %>%
  arrange(Scenario, Habitat_from, Species)
write_csv(Species, Out_species_csv)
Locations <- map(Scenario_results, "locations") %>%
  list_rbind() %>%
  mutate(across(c(starts_with("dRich"), Cover_pct, all_of(Ssp_col)), \(x) round(x, 3)), Tot_prec = round(Tot_prec, 1)) %>%
  arrange(Scenario, Site)
write_csv(Locations, Out_locations_csv)
Landscape <- map(Scenario_results, "landscape") %>%
  list_rbind() %>%
  mutate(across(c(Rich_landscape_before, starts_with("dRich")), \(x) round(x, 3)), P_gt0 = round(P_gt0, 3))
write_csv(Landscape, Out_landscape_csv)
dir.create(Fig_dir, showWarnings = FALSE, recursive = TRUE)

## Change in expected richness by scenario and original habitat: posterior median and CrI of the mean over changed location-periods, with each changed location's posterior median behind it (buffer scenarios with canopy cover changed)
Main_scenarios <- c("disc_to_ssp", "disc_to_pasture", "buffer_to_ssp", "buffer_to_pasture")
Plot_df <- filter(Summary, Habitat_from != "All", Scenario %in% Main_scenarios) %>% mutate(Scenario_key = Scenario, Scenario = factor(Scenario_labels[Scenario], levels = Scenario_labels[Main_scenarios]))
Loc_df <- filter(Locations, Scenario %in% Main_scenarios) %>% mutate(Scenario = factor(Scenario_labels[Scenario], levels = Scenario_labels[Main_scenarios]))
Rich_label_y <- max(Loc_df$dRich_median, Plot_df$dRich_high) + 2
## Two panels, one per conversion extent (a: point count circle, b: all pasture or silvopasture within 300 m), on a shared y-axis with each habitat's median change written in one row above the data, so the habitats' changes within a scenario can be read off and added together
## A second version gives, in panel b, the overall change in landscape richness ignoring the sampling design (union of the three habitats' mean communities) in a row below the data; both panels of a version share its y-axis
Overall_label_y <- min(Loc_df$dRich_median, Plot_df$dRich_low) - 2.5
Rich_limits <- range(Loc_df$dRich_median, Plot_df$dRich_low, Rich_label_y)
Rich_limits_overall <- range(Rich_limits, Overall_label_y)
signed <- \(x) str_replace(sprintf("%+.1f", x), "-", "−")
Direction_labels <- c(disc_to_ssp = "Pasture to\nsilvopasture", disc_to_pasture = "Silvopasture\nto pasture", buffer_to_ssp = "All pasture to\nsilvopasture", buffer_to_pasture = "All silvopasture\nto pasture")
conversion_panel <- function(scenarios, title, y_title = "Change in expected species richness", overall = FALSE, y_limits = Rich_limits) {
  as_direction <- \(df) filter(df, Scenario_key %in% scenarios) %>% mutate(Scenario = factor(Direction_labels[Scenario_key], levels = Direction_labels[scenarios]))
  ggplot(as_direction(Plot_df), aes(x = Scenario, y = dRich_median, color = Habitat_from)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(data = as_direction(mutate(Loc_df, Scenario_key = Main_scenarios[match(Scenario, Scenario_labels[Main_scenarios])])), aes(y = dRich_median), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 20261107), alpha = 0.3, size = 1) +
    geom_pointrange(aes(ymin = dRich_low, ymax = dRich_high), position = position_dodge(width = 0.8), linewidth = 1, size = 0.6) +
    geom_text(aes(y = Rich_label_y, label = signed(dRich_median)), position = position_dodge(width = 0.8), size = 4.3, fontface = "bold", show.legend = FALSE) +
    {if (overall) geom_text(data = as_direction(mutate(Landscape, Scenario_key = Scenario)), aes(x = Scenario, y = Overall_label_y, label = paste0("Overall: ", signed(dRich_landscape_median), " (", signed(dRich_landscape_low), ", ", signed(dRich_landscape_high), ")")), inherit.aes = FALSE, size = 4.3, fontface = "bold", color = "grey20")} +
    scale_color_manual(values = level_colors, limits = names(level_colors), name = "Habitat before") +
    scale_y_continuous(limits = y_limits) +
    labs(x = NULL, y = y_title, title = title) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"))
}
Conversion_plots <- map(c(plain = FALSE, overall = TRUE), \(overall) {
  y_limits <- if (overall) Rich_limits_overall else Rich_limits
  patchwork::wrap_plots(conversion_panel(Main_scenarios[1:2], "Point count circle (50 m)", y_limits = y_limits) + theme(legend.position = "none"), conversion_panel(Main_scenarios[3:4], "Within 300 m", y_title = NULL, overall = overall, y_limits = y_limits), widths = c(2, 3), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
})
ggsave(file.path(Fig_dir, "Habitat_conversion_richness.png"), Conversion_plots$plain, width = 11, height = 5.5, bg = "white")
print(Conversion_plots$plain)
ggsave(file.path(Fig_dir, "Habitat_conversion_richness_overall.png"), Conversion_plots$overall, width = 11, height = 5.5, bg = "white")
print(Conversion_plots$overall)

## Expected richness before (open) and after (filled, CrI) conversion by scenario and original habitat, against the mean expected richness of the real locations of each habitat (dashed), labelled with the share of the gap to the new habitat closed (buffer scenarios with canopy cover changed)
Habitat_rich_df <- tibble(Habitat = names(Rich_by_habitat), Rich = map_dbl(Rich_by_habitat, median))
Before_after_plot <- ggplot(Plot_df, aes(x = Habitat_from, color = Habitat_from)) +
  geom_hline(data = Habitat_rich_df, aes(yintercept = Rich, color = Habitat), linetype = "dashed") +
  geom_segment(aes(xend = Habitat_from, y = Rich_obs, yend = Rich_after_median), arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(y = Rich_obs), shape = 1, size = 2.5) +
  geom_pointrange(aes(y = Rich_after_median, ymin = Rich_after_low, ymax = Rich_after_high), linewidth = 1, size = 0.5) +
  geom_text(data = filter(Plot_df, !is.na(Gap_closed_median)), aes(y = Rich_after_high, label = sprintf("%.0f%%", 100 * Gap_closed_median)), vjust = -0.6, size = 3.5, show.legend = FALSE) +
  facet_wrap(vars(Scenario), nrow = 1) +
  scale_color_manual(values = level_colors, name = "Habitat") +
  labs(x = "Habitat before", y = "Expected species richness") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(Fig_dir, "Habitat_conversion_richness_before_after.png"), Before_after_plot, width = 11, height = 5, bg = "white")
print(Before_after_plot)

## Species: change in occupancy at the converted points of the two 50 m disc scenarios, ordered by the pasture-to-silvopasture change
Sp_df <- filter(Species, Scenario %in% c("disc_to_ssp", "disc_to_pasture")) %>%
  mutate(Scenario = factor(Scenario_labels[Scenario], levels = Scenario_labels[c("disc_to_ssp", "disc_to_pasture")]), Species = str_replace(Species, "_", " "))
Sp_order <- filter(Sp_df, Scenario == Scenario_labels[["disc_to_ssp"]]) %>% arrange(dPsi_median) %>% pull(Species)
Species_plot <- ggplot(mutate(Sp_df, Species = factor(Species, levels = Sp_order)), aes(x = dPsi_median, y = Species, color = Habitat_from)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = dPsi_low, xmax = dPsi_high), size = 0.15, linewidth = 0.4) +
  facet_wrap(vars(Scenario)) +
  scale_color_manual(values = level_colors, name = "Habitat before") +
  labs(x = "Change in occupancy probability", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic", size = 6), legend.position = "bottom")
ggsave(file.path(Fig_dir, "Habitat_conversion_species.png"), Species_plot, width = 8, height = 12, bg = "white")
print(Species_plot)

## Context dependence: change at each converted location in the 50 m disc scenarios against its canopy cover within 1.5 km, silvopasture within 300 m and rainfall, with a straight line through the location medians of each direction
Context_df <- filter(Locations, Scenario %in% c("disc_to_ssp", "disc_to_pasture")) %>%
  pivot_longer(c(Cover_pct, all_of(Ssp_col), Tot_prec), names_to = "Covariate", values_to = "Value") %>%
  mutate(Covariate = factor(c(Cover_pct = "Canopy cover within 1.5 km (%)", set_names(Ssp_label, Ssp_col), Tot_prec = "Annual precipitation (mm)")[Covariate], levels = c("Canopy cover within 1.5 km (%)", Ssp_label, "Annual precipitation (mm)")))
Context_plot <- ggplot(Context_df, aes(x = Value, y = dRich_median, color = Habitat_from)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = dRich_low, ymax = dRich_high), size = 0.2, linewidth = 0.3, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9) +
  facet_wrap(vars(Covariate), scales = "free_x", strip.position = "bottom") +
  scale_color_manual(values = level_colors, labels = c(Pasture = "Pasture to silvopasture", Silvopasture = "Silvopasture to pasture"), name = "Point count circle\n(50 m) converted") +
  labs(x = NULL, y = "Change in expected species richness") +
  theme_bw() +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = "bottom")
ggsave(file.path(Fig_dir, "Habitat_conversion_context.png"), Context_plot, width = 10, height = 4.5, bg = "white")
print(Context_plot)

# Composition ordination ----
## NMDS of the real location-periods on the posterior-mean Bray-Curtis dissimilarity of their fitted occupancy (each ordination draw's matrix averaged element-wise); converted location-periods are added from their mean dissimilarity to every real location-period over the same draws; MDSaddpoints needs exactly the same index for both, so both come from one Bray-Curtis formula (designdist / designdist2, minimum terms; designdist2(x, y) returns rows of y against columns of x)
Bray_formula <- "(A+B-2*J)/(A+B)"
Dist_real <- map(Ord_draws, \(d) as.matrix(vegan::designdist(t(Psi_obs[d, , ]), Bray_formula, terms = "minimum"))) %>%
  reduce(`+`) / Ord_n_draws
set.seed(Ord_seed)
Nmds_real <- vegan::metaMDS(as.dist(Dist_real), k = Ord_k, try = Ord_try, trymax = Ord_try, trace = 0)
Real_scores <- tibble(Cell = seq_len(nrow(Cells)), Habitat = unname(level_labels[Cells$Habitat]), NMDS1 = Nmds_real$points[, 1], NMDS2 = Nmds_real$points[, 2])
Added_scores <- map(Main_scenarios, \(scenario) {
  res <- Scenario_results[[scenario]]
  dis <- map(seq_len(Ord_n_draws), \(k) vegan::designdist2(t(Psi_obs[Ord_draws[k], , ]), t(res$psi_new[k, , ]), Bray_formula, terms = "minimum")) %>%
    reduce(`+`) / Ord_n_draws
  added <- vegan::MDSaddpoints(Nmds_real, unclass(dis))
  tibble(Scenario = scenario, Cell = res$changed, NMDS1_new = added$points[, 1], NMDS2_new = added$points[, 2])
}) %>%
  list_rbind() %>%
  left_join(Real_scores, by = "Cell", relationship = "many-to-one") %>%
  rename(Habitat_from = Habitat) %>%
  mutate(Scenario = factor(Scenario_labels[Scenario], levels = Scenario_labels[Main_scenarios]))
## 1-SD dispersion ellipse of each habitat's real location-periods on the first two axes (as vegan::ordiellipse(kind = "sd"))
Ellipses <- Real_scores %>%
  group_by(Habitat) %>%
  group_modify(\(df, key) {
    xy <- cbind(df$NMDS1, df$NMDS2)
    circle <- cbind(cos(seq(0, 2 * pi, length.out = 100)), sin(seq(0, 2 * pi, length.out = 100)))
    pts <- sweep(circle %*% chol(cov(xy)), 2, colMeans(xy), "+")
    tibble(NMDS1 = pts[, 1], NMDS2 = pts[, 2])
  }) %>%
  ungroup()
Added_centroids <- summarise(Added_scores, across(c(NMDS1, NMDS2, NMDS1_new, NMDS2_new), mean), .by = c(Scenario, Habitat_from))

## Real location-periods with their habitat's 1-SD ellipse, each converted location-period's move (thin arrows) and the move of each original habitat's mean position (bold arrows), by scenario (buffer scenarios with canopy cover changed)
Nmds_plot <- ggplot(mapping = aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = Real_scores, aes(color = Habitat), size = 0.6, alpha = 0.25) +
  geom_path(data = Ellipses, aes(color = Habitat), linewidth = 0.7) +
  geom_segment(data = Added_scores, aes(xend = NMDS1_new, yend = NMDS2_new, color = Habitat_from), linewidth = 0.25, alpha = 0.35, arrow = arrow(length = unit(0.08, "cm"))) +
  geom_segment(data = Added_centroids, aes(xend = NMDS1_new, yend = NMDS2_new, color = Habitat_from), linewidth = 1.3, arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  facet_wrap(vars(Scenario), nrow = 2) +
  scale_color_manual(values = level_colors, name = "Habitat") +
  coord_equal() +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")
ggsave(file.path(Fig_dir, "Habitat_conversion_nmds.png"), Nmds_plot, width = 10, height = 6.4, bg = "white")
print(Nmds_plot)

## What the first two axes track: Spearman correlation of the real location-periods' scores with their covariates (scaled model columns, or percent of the 300 m buffer) and expected richness, over all location-periods and within each habitat; and the species with the lowest and highest weighted-average scores on axis 2 (weights = posterior-mean occupancy)
Psi_mean <- t(colMeans(Psi_obs))
Axis_covs <- tibble(Habitat = Real_scores$Habitat, NMDS1 = Real_scores$NMDS1, NMDS2 = Real_scores$NMDS2, Richness = rowSums(Psi_mean), !!Ssp_col := Cells$Ssp_fit, Forest_300 = Cells$forest, Pasture_300 = Cells$intpast,
  Cover_1500 = Cells$Cover_pct, Height_1500 = X_obs[, "Height_h1500"], Tot_prec = X_obs[, "Tot_prec"], Dist_open = X_obs[, "Dist_open"], Dist_inside = X_obs[, "Dist_inside"], Year = X_obs[, "Year_cont"])
Axis_cor <- map(c(list(All = Axis_covs), split(Axis_covs, Axis_covs$Habitat)), \(df) {
  covs <- select(df, -Habitat, -NMDS1, -NMDS2) %>% select(where(\(x) sd(x) > 0))
  tibble(Covariate = names(covs), NMDS1 = round(cor(df$NMDS1, covs, method = "spearman")[1, ], 2), NMDS2 = round(cor(df$NMDS2, covs, method = "spearman")[1, ], 2))
}) %>% list_rbind(names_to = "Location_periods")
Species_scores <- vegan::wascores(Nmds_real$points[, 1:2], Psi_mean) %>%
  as_tibble(rownames = "Species") %>%
  mutate(Species = sp_names, Mean_psi = colMeans(Psi_mean)) %>%
  arrange(MDS2)

## Do the habitats' communities converge? Bray-Curtis between the mean communities of each pair of original habitats (all their location-periods, converted ones at their counterfactual occupancy), before and after each scenario, at the ordination draws
Habitat_pairs <- combn(names(Habitat_cells), 2, simplify = FALSE)
Convergence <- map(Main_scenarios, \(scenario) {
  res <- Scenario_results[[scenario]]
  psi_after <- Psi_obs[Ord_draws, , , drop = FALSE]
  psi_after[, , res$changed] <- res$psi_new
  map(Habitat_pairs, \(pair) {
    before <- bray_curtis_rows(rowMeans(Psi_obs[Ord_draws, , Habitat_cells[[pair[1]]], drop = FALSE], dims = 2), rowMeans(Psi_obs[Ord_draws, , Habitat_cells[[pair[2]]], drop = FALSE], dims = 2))
    after <- bray_curtis_rows(rowMeans(psi_after[, , Habitat_cells[[pair[1]]], drop = FALSE], dims = 2), rowMeans(psi_after[, , Habitat_cells[[pair[2]]], drop = FALSE], dims = 2))
    q <- quantile(after - before, probs)
    tibble(Scenario = scenario, Pair = paste(pair, collapse = " - "), BC_before = round(median(before), 3), BC_after = round(median(after), 3), dBC_low = round(q[[1]], 3), dBC_median = round(q[[2]], 3), dBC_high = round(q[[3]], 3), P_closer = round(mean(after < before), 3))
  }) %>% list_rbind()
}) %>% list_rbind()

## Console report
cat("Composition NMDS of the real location-periods: stress", round(Nmds_real$stress, 3), "in", Ord_k, "dimensions\n")
cat("Spearman correlation of the NMDS axes with covariates:\n"); print(pivot_wider(Axis_cor, names_from = Location_periods, values_from = c(NMDS1, NMDS2)), width = Inf)
cat("Species at the low and high ends of NMDS axis 2:\n"); print(bind_rows(head(Species_scores, 8), tail(Species_scores, 8)), n = 16)
cat("Bray-Curtis between the mean communities of the original habitats, before and after each scenario:\n"); print(Convergence, n = Inf, width = Inf)
cat("Disc share of the buffer:", round(Disc_share_pct, 2), "percentage points\n")
cat("Canopy cover in the 50 m disc (%):", paste(Disc_by_habitat$Habitat, round(Disc_by_habitat$Canopy_disc, 3), "( n =", Disc_by_habitat$N, ")", collapse = "; "), "| gain on conversion", round(Canopy_gain, 3), "\n")
cat("Ring canopy cover change (percentage points), all pasture to silvopasture:", round(quantile(Cells$dCover_to_ssp, c(0, 0.5, 1)), 2), "| all silvopasture to pasture:", round(quantile(Cells$dCover_to_pasture, c(0, 0.5, 1)), 2), "\n")
cat("Mean land cover of the 50 m disc (%):\n"); print(summarise(Cells, across(c(forest_50, intpast_50, ssp_50), \(x) round(mean(x), 1)), .by = Habitat))
print(Ssp_max)
cat("Canopy density of farm silvopasture (canopy-weighted / area-based), over", nrow(Farm_density), "farms:", round(quantile(Farm_density$Density, c(0, 0.25, 0.5, 0.75, 1)), 2), "\n")
print(Summary, width = Inf)

# Scenario schematic ----
## Diagram of the two conversion extents (pasture to silvopasture shown; the reverse scenarios mirror it) on an invented landscape, not data: one survey point per habitat with its 50 m point count circle (solid) and 300 m buffer (dashed), drawn to scale; the baseline landscape on the left, arrows to the two alternative simulations stacked on the right
set.seed(20261107)
Schematic_pts <- tibble(Habitat = c("Pasture", "Silvopasture", "Forest"), x = c(330, 800, 1260), y = 430)
Pts_sfc <- sf::st_sfc(map2(Schematic_pts$x, Schematic_pts$y, \(x, y) sf::st_point(c(x, y))))
Landscape <- sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1600, ymax = 860)))
## Irregular closed outline for a forest fragment, centred at (x, y) with radius about r (the last vertex repeats the first so the ring closes exactly)
blob <- function(x, y, r) {
  angle <- seq(0, 2 * pi, length.out = 40)[-40]
  radius <- r * (1 + 0.18 * sin(3 * angle + 1) + 0.12 * cos(5 * angle))
  xy <- cbind(x + radius * cos(angle), y + radius * sin(angle))
  sf::st_polygon(list(rbind(xy, xy[1, ])))
}
## Forest: a riparian strip and three fragments, two reaching into the pasture and forest points' buffers
Forest_poly <- c(sf::st_buffer(sf::st_sfc(sf::st_linestring(cbind(c(1230, 1310, 1230, 1300), c(0, 290, 580, 860)))), 85), sf::st_sfc(blob(150, 175, 105), blob(560, 755, 80), blob(1475, 170, 85))) %>% sf::st_union() %>% sf::st_intersection(Landscape)
## Silvopasture: a dispersed-tree pasture around the silvopasture point and live fences (tree rows) along field boundaries, some crossing the buffers, each ending at the forest edge
Fence_lines <- sf::st_sfc(sf::st_linestring(cbind(c(0, 640), c(610, 610))), sf::st_linestring(cbind(c(560, 560), c(0, 700))), sf::st_linestring(cbind(c(980, 1600), c(150, 150))), sf::st_linestring(cbind(c(1040, 1040), c(640, 860))), sf::st_linestring(cbind(c(1350, 1600), c(560, 560)))) %>%
  sf::st_difference(sf::st_buffer(Forest_poly, 15)) %>%
  sf::st_cast("MULTILINESTRING") %>%
  sf::st_cast("LINESTRING")
Fence_lines <- Fence_lines[as.numeric(sf::st_length(Fence_lines)) > 40]
Ssp_block <- sf::st_buffer(sf::st_sfc(sf::st_polygon(list(cbind(c(660, 960, 990, 700, 660), c(250, 230, 620, 650, 250))))), 25)
Ssp_poly <- c(Ssp_block, sf::st_buffer(Fence_lines, 12, endCapStyle = "FLAT")) %>% sf::st_intersection(Landscape)
## Isolated trees in the open pasture, drawn on top of the land cover as small silvopasture-coloured dots (the live fences and dispersed-tree pasture are shown as silvopasture area only)
Trees <- sf::st_sample(sf::st_difference(Landscape, sf::st_union(c(Forest_poly, sf::st_buffer(sf::st_union(Ssp_poly), 30), sf::st_buffer(Pts_sfc, 70)))), 22, type = "random")
Circles_50 <- sf::st_buffer(Pts_sfc, 50)
Buffers_300 <- sf::st_union(sf::st_buffer(Pts_sfc, 300)) %>% sf::st_intersection(Landscape)
## Panel layout: baseline landscape on the left, centred on a 2 x 2 grid of simulations (columns: conversion direction; rows: conversion extent)
Panel_w <- 1600; Panel_h <- 860; Panel_gap <- 260; Grid_x <- 2700
Schematic_layout <- tibble(Panel = c("Real", "Ssp_circle", "Ssp_300", "Pasture_circle", "Pasture_300"),
  dx = c(0, Grid_x, Grid_x, Grid_x + Panel_w + Panel_gap, Grid_x + Panel_w + Panel_gap), dy = c((Panel_h + Panel_gap) / 2, Panel_h + Panel_gap, 0, Panel_h + Panel_gap, 0))
shift <- \(geom, panel) geom + c(Schematic_layout$dx[Schematic_layout$Panel == panel], Schematic_layout$dy[Schematic_layout$Panel == panel])
## Land cover of one panel from its silvopasture area: forest is never converted, pasture is the rest of the landscape
schematic_cover <- function(ssp_area, panel) {
  ssp_area <- sf::st_difference(sf::st_union(ssp_area), Forest_poly)
  sf::st_sf(Cover = factor(c("Pasture", "Silvopasture", "Forest"), levels = names(level_colors)), geometry = shift(c(sf::st_difference(Landscape, sf::st_union(ssp_area, Forest_poly)), ssp_area, Forest_poly), panel))
}
## Silvopasture area and remaining isolated trees in each panel: pasture to silvopasture adds the pasture point's circle or every point's 300 m buffer (isolated trees inside it are absorbed into the silvopasture); silvopasture to pasture removes the silvopasture point's circle or all silvopasture (and isolated trees) within every 300 m buffer
Ssp_union <- sf::st_union(Ssp_poly)
Panel_ssp <- list(Real = Ssp_union, Ssp_circle = sf::st_union(Ssp_union, Circles_50[1]), Ssp_300 = sf::st_union(Ssp_union, Buffers_300), Pasture_circle = sf::st_difference(Ssp_union, Circles_50[2]), Pasture_300 = sf::st_difference(Ssp_union, Buffers_300))
Panel_trees <- list(Real = Trees, Ssp_circle = Trees, Ssp_300 = Trees[!lengths(sf::st_intersects(Trees, Buffers_300))], Pasture_circle = Trees[!lengths(sf::st_intersects(Trees, Circles_50[2]))], Pasture_300 = Trees[!lengths(sf::st_intersects(Trees, Buffers_300))])
Schematic_cover <- imap(Panel_ssp, \(ssp_area, panel) schematic_cover(ssp_area, panel)) %>% bind_rows()
Schematic_trees <- imap(Panel_trees, \(trees, panel) as_tibble(sf::st_coordinates(shift(trees, panel)))) %>% list_rbind()
in_panels <- \(geom) map(Schematic_layout$Panel, \(panel) shift(geom, panel)) %>% do.call(c, .)
Schematic_rings <- sf::st_sf(Ring = rep(c("50 m", "300 m"), each = 3 * nrow(Schematic_layout)), geometry = c(in_panels(Circles_50), in_panels(sf::st_buffer(Pts_sfc, 300))))
Schematic_points <- sf::st_coordinates(in_panels(Pts_sfc)) %>% as_tibble()
Schematic_labels <- mutate(Schematic_pts, x = x + Schematic_layout$dx[1], y = y + Schematic_layout$dy[1] + 90, Label = Habitat)
## Titles: the baseline landscape, the two conversion directions above the grid columns, and the two conversion extents along the arrows to the grid rows
Schematic_titles <- tibble(x = c(Panel_w / 2, Grid_x + Panel_w / 2, Grid_x + 1.5 * Panel_w + Panel_gap), y = c(Schematic_layout$dy[1] + Panel_h + 70, rep(2 * Panel_h + Panel_gap + 70, 2)), Label = c("Baseline landscape", "Pasture → silvopasture", "Silvopasture → pasture"))
Schematic_arrows <- tibble(x = Panel_w + 60, y = Schematic_layout$dy[1] + Panel_h / 2, xend = Grid_x - 60, yend = c(Panel_h + Panel_gap, 0) + Panel_h / 2, Label = c("Point count\ncircle (50 m)", "All within\n300 m"))
Schematic_plot <- ggplot() +
  geom_sf(data = Schematic_cover, aes(fill = Cover), color = NA, alpha = 0.6) +
  geom_point(data = Schematic_trees, aes(X, Y), shape = 21, size = 1.6, fill = level_colors[["Silvopasture"]], color = "white", stroke = 0.2) +
  geom_sf(data = Schematic_rings, aes(linetype = Ring), fill = NA, color = "grey20", linewidth = 0.5, key_glyph = "path") +
  geom_point(data = Schematic_points, aes(X, Y), size = 1.6) +
  geom_label(data = Schematic_labels, aes(x, y, label = Label), size = 3.6, label.size = 0, fill = alpha("white", 0.8)) +
  geom_text(data = Schematic_titles, aes(x, y, label = Label), size = 5.5, fontface = "bold", vjust = 0) +
  geom_segment(data = Schematic_arrows, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(0.18, "inches"), type = "closed"), linewidth = 1, color = "grey25") +
  geom_text(data = Schematic_arrows, aes((x + xend) / 2 - 100, (y + yend) / 2 + c(1, -1) * 300, label = Label), size = 4.6, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = level_colors, name = "Land cover") +
  scale_linetype_manual(values = c(`50 m` = "solid", `300 m` = "dashed"), breaks = c("50 m", "300 m"), name = "Around each point") +
  coord_sf(datum = NA, clip = "off", xlim = c(-20, Grid_x + 2 * Panel_w + Panel_gap + 20), ylim = c(-20, 2 * Panel_h + Panel_gap + 200), expand = FALSE) +
  theme_void(base_size = 14)
ggsave(file.path(Fig_dir, "Habitat_conversion_schematic.png"), Schematic_plot, width = 16, height = 5.6, bg = "white")
print(Schematic_plot)
