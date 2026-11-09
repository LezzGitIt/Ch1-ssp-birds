## PhD birds in silvopastoral landscapes ##
# Scale of effect of landscape canopy cover on occupancy: extracts woody canopy cover (trees / shrubs > 2 m, Colombia Woody Vegetation Structure and Change, ~27 m) in rings around each point count over a wide range of radii, matched to the survey year, then checks how the covariate behaves across radii (variation, correlation among radii, collinearity with the other occupancy covariates) before any model is fitted.

### Rings, not discs: each ring excludes the inner Ring_inner_m around the point, so the landscape measure does not duplicate the point's own Habitat (as in Karp et al. 2019).
### Only the surveyed location x period rows of the stacked data list (the cells that inform the multi-season model) enter the diagnostics.
### Scale-of-effect fits (last section, run only when named on the command line): the primary model with one ring's canopy cover added, one fit per radius, plus the primary model refitted unchanged as the baseline.

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

## Outputs
Rings_csv <- "Derived/Excels/Canopy_rings_by_location.csv"
Fig_dir <- "Figures/Scale_of_effect"

## Fits to run, from the command line: 'baseline', 'r<radius>' (e.g. r300), or 'smoke' (10 batches, 1 chain, with the 1 km ring; checks the terms and MCMC lengths against the primary fit's summary); none when the script is sourced
run_sections <- commandArgs(trailingOnly = TRUE)

## Primary model (Multi_spp_100 row 90), copied from spOccupancy_temporal.R's configuration and its 2-factor section; every fit here is this model plus at most one ring term, all with the same seed
tmsOcc_l_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_98_hab_fix.rds"
Primary_summary_path <- "Derived/Rdata/Model_output/Temporal/Occ_mod_98spp_prec_nofarm_long_f2_2026-11-05_summary.rds"
Primary_occ_formula <- ~ Habitat * Tot_prec + Dist_forest + ssp + Year1
Primary_fit_settings <- list(
  det.formula = ~ Forest_bin + Pc_start + I(Pc_start^2) + Prob_breed + Sampling_day_log + (1 | Collector_team_num),
  priors = list(phi.unif = list(3 / 25e3, 3 / 23e3)),
  n.neighbors = 5, n.factors = 2, batch.length = 25, ar1 = FALSE,
  mcmc = list(n.batch = 2400, n.burn = 30000, n.thin = 75, n.chains = 3),
  out_dir = "Derived/Rdata/Model_output/Temporal/Scale_of_effect",
  seed = 20261107
)

# Load data ----
## Surveyed rows of the model, with the survey year taken as in the wrangling script (latest year within a location x year-group x season) and snapped to the nearest raster year (2025-26 surveys use 2024)
Event_years <- read_csv(paste0(Wrangling_excels, "Event_covs.csv"), show_col_types = FALSE) %>%
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
        terra::crop(terra::rast(file.path(National_dir, sprintf("Colombia_WVCC_%d.tif", grp$Raster_year[1]))), terra::buffer(pts, max(Radii_m) + 500), snap = "out")
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

  Canopy_rings <- Disc_sums %>%
    filter(Radius_m != Ring_inner_m) %>%
    left_join(Disc_sums %>% filter(Radius_m == Ring_inner_m) %>% select(Id_survey_no_dc, Raster_year, sum_inner = sum, count_inner = count), by = c("Id_survey_no_dc", "Raster_year")) %>%
    mutate(Canopy_cover_ring = round((sum - sum_inner) / (count - count_inner), 3)) %>%
    select(Id_survey_no_dc, Raster_year, Radius_m, Canopy_cover_ring) %>%
    arrange(Id_survey_no_dc, Raster_year, Radius_m)
  write_csv(Canopy_rings, Rings_csv)
}
stopifnot(!anyNA(Canopy_rings$Canopy_cover_ring))

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

# Report ----
cat(nrow(Rows), "surveyed rows,", nrow(Loc_years), "location x raster-year combinations,", n_distinct(Rows$Id_scr), "farms\n")
print(Scale_summary %>% mutate(across(-Radius_m, \(x) round(x, 3))), width = Inf)

# Scale-of-effect fits ----
### Each fit is the primary model plus one ring's canopy cover as an additive term ('r<radius>'), or the primary model unchanged ('baseline'); radii are then compared by WAIC relative to the baseline and by the canopy coefficient. The baseline and the unseeded row 90 fit are two runs of one model, so their WAIC difference shows how much WAIC moves between runs.
### Run headless, about 3 h per fit, two processes at a time: caffeinate -is Rscript Scripts/Scale_of_effect.R r150 r300 > <log> 2>&1
Soe_fit_tags <- intersect(run_sections, c("smoke", "baseline", paste0("r", Radii_m)))
if (length(Soe_fit_tags) > 0) {
  primary_summary <- readRDS(Primary_summary_path)
  mcmc <- Primary_fit_settings$mcmc
  stopifnot(all(c(mcmc$n.batch * Primary_fit_settings$batch.length, mcmc$n.burn, mcmc$n.thin, mcmc$n.chains) == primary_summary$mcmc[c("n.samples", "n.burn", "n.thin", "n.chains")]))

  ## Ring canopy cover as one site x period matrix per radius, in the multi-season list's site and period order: standardised over the surveyed rows (as the wrangling script scales every covariate), unsurveyed cells filled from the location's nearest surveyed period (as the list's other time-varying covariates)
  tmsOcc_l <- readRDS(tmsOcc_l_path)
  Site_levels <- dimnames(tmsOcc_l$y)[[2]]
  Time_levels <- dimnames(tmsOcc_l$y)[[3]]
  Surveyed <- apply(!is.na(tmsOcc_l$y[1, , , , drop = FALSE]), c(2, 3), any)
  Ring_covs <- Rows_rings %>%
    mutate(Canopy_ring = as.numeric(scale(Canopy_cover_ring)), .by = Radius_m) %>%
    split(.$Radius_m) %>%
    map(\(d) {
      m <- matrix(NA_real_, length(Site_levels), length(Time_levels), dimnames = list(Site_levels, Time_levels))
      m[cbind(match(d$Site_code, Site_levels), match(paste(d$Year_grp, d$Season, sep = "."), Time_levels))] <- d$Canopy_ring
      stopifnot(identical(unname(!is.na(m)), unname(matrix(Surveyed, nrow(m)))))
      fill_nearest_period(m)
    }) %>%
    set_names(\(radius) paste0("Canopy_r", radius))
  tmsOcc_l$occ.covs <- c(tmsOcc_l$occ.covs, Ring_covs)

  walk(Soe_fit_tags, \(fit_tag) {
    ring_term <- switch(fit_tag, baseline = NULL, smoke = "Canopy_r1000", paste0("Canopy_", fit_tag))
    occ_formula <- if (is.null(ring_term)) Primary_occ_formula else update(Primary_occ_formula, as.formula(paste("~ . +", ring_term)))
    mod <- fit_temporal(tmsOcc_l, occ_formula, paste0("soe_", fit_tag), settings = Primary_fit_settings, smoke = fit_tag == "smoke")
    # Same occupancy terms as the primary fit plus the ring term, and the same detection terms
    stopifnot(setequal(colnames(mod$beta.comm.samples), c(primary_summary$occ_terms, ring_term)), setequal(colnames(mod$alpha.comm.samples), primary_summary$det_terms))
    if (fit_tag == "smoke") cat("Smoke fit OK:", deparse(occ_formula), "\n")
  })
}
