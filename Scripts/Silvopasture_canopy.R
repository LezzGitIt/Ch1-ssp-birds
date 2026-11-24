## PhD birds in silvopastoral landscapes ##
# Canopy-weighted silvopasture cover within 300 m of each point count, and the number of digitized isolated trees within 300 m, from the digitized Meta land-cover polygons; the area-based silvopasture cover (ssp) counts a sparse dispersed-tree stand like a dense one, which overstates dispersed-tree canopy about threefold.

### Canopy weights by silvopasture type (user, 2026-11-23): dispersed-tree stands at their imagery-based tree cover (trees_perc), live fences, intensive silvopasture and isolated-tree crowns (points buffered 4 m) at 100%, forage banks excluded (0%).
### Each land-cover period's polygons are matched to the surveys as Farm_diversity's landscape-metrics script matches its ssp (the 'middle' layer by default, overwritten by the 'past' and 'ubc' layers in their data year, gaps filled from the same point's other years), and the area-based ssp is carried through the same steps as a check that the mapping is reproduced.

# Load libraries ----
library(sf)
library(tidyverse)
library(readxl)
library(cowplot)
ggplot2::theme_set(cowplot::theme_cowplot())

# Configuration ----
Lc_dir <- "../Farm_diversity/Derived/Geospatial/shp/R_processed"
Lsm_dir <- "../Farm_diversity/Derived/Excels/Lsm"
Wrangling_dir <- "../Ssp-bird-data-wrangling/Derived/"
## Point-digitized isolated trees (one layer for all periods and regions, no date attribute)
Trees_path <- "../../Mentorship/Digitization_Mathilde/Final_docs/final_shp_for_natalia/current_docs/individual_trees_final"
Lc_periods <- c("middle", "past", "ubc")
Buffer_m <- 300
Crown_radius_m <- 4
## A dispersed-tree polygon counts as an isolated-tree crown when it has no trees_perc (crowns, crown fragments and pairs of crowns merged during snapping, 72-110 m2, 50-99% inside the crown buffers) or when at least this share of it lies within the crown radius of a digitized tree (crowns in the 'ubc' layer inherit the trees_perc of the polygon they were cut from)
Crown_min_share <- 0.99
Canopy_weights <- c(livefence = 1, silvopast = 1, forbanks = 0)
Out_csv <- "Derived/Excels/Ssp_canopy_by_survey.csv"
Trees_csv <- "Derived/Excels/Isolated_trees_by_location.csv"
Fig_dir <- "Figures/Silvopasture_canopy"
dir.create(Fig_dir, showWarnings = FALSE)

# Load data ----
Polys <- map(Lc_periods, \(p) st_read(file.path(Lc_dir, paste0("Snapped_", p, ".gpkg")), quiet = TRUE) |> select(Id_survey, Department, lc_typ, lc_typ2, trees_perc) |> mutate(Lc_period = p)) |>
  bind_rows() |>
  filter(Department == "Meta") |>
  mutate(Poly_id = row_number(), Area_m2 = as.numeric(st_area(geom)))
Trees <- st_read(Trees_path, quiet = TRUE) |> st_transform(st_crs(Polys))
Pc_locs <- st_read(paste0(Wrangling_dir, "Geospatial/shp/Pc_locs_dc.gpkg"), quiet = TRUE) |>
  filter(Department == "Meta") |>
  st_transform(st_crs(Polys))

# Isolated-tree crowns ----
## Share of each dispersed-tree polygon inside the union of the 4 m crowns of the digitized trees
Crowns <- st_union(st_buffer(Trees, Crown_radius_m + 0.05))
Crown_share <- Polys |>
  filter(lc_typ == "distrees") |>
  st_make_valid() |>
  st_intersection(Crowns) |>
  mutate(Crown_m2 = as.numeric(st_area(geom))) |>
  st_drop_geometry() |>
  summarise(Crown_m2 = sum(Crown_m2), .by = Poly_id)

# Canopy-weighted ssp ----
Poly_weights <- Polys |>
  st_drop_geometry() |>
  left_join(Crown_share, by = "Poly_id") |>
  mutate(Crown = lc_typ == "distrees" & (is.na(trees_perc) | coalesce(Crown_m2 / Area_m2, 0) >= Crown_min_share),
         Canopy_weight = case_when(lc_typ2 != "ssp" ~ 0, Crown ~ 1, lc_typ == "distrees" ~ trees_perc / 100, .default = Canopy_weights[lc_typ]))
cat("Dispersed-tree polygons by crown status and trees_perc:\n"); print(count(filter(Poly_weights, lc_typ == "distrees"), Lc_period, Crown, Missing_trees_perc = is.na(trees_perc)))
stopifnot(!anyNA(Poly_weights$Canopy_weight))

## Shares of all digitized area within the 300 m buffer (%), as the landscape metrics' ssp is a share of the digitized landscape
Ssp_period <- Poly_weights |>
  summarise(Ssp_area_poly = 100 * sum(Area_m2 * (lc_typ2 == "ssp")) / sum(Area_m2),
            Ssp_cw = 100 * sum(Area_m2 * Canopy_weight) / sum(Area_m2),
            Crown_pct = 100 * sum(Area_m2 * Crown) / sum(Area_m2),
            Stand_pct = 100 * sum(Area_m2 * (lc_typ == "distrees" & !Crown)) / sum(Area_m2),
            N_small_stands = sum(lc_typ == "distrees" & !Crown & Area_m2 < 500), .by = c(Lc_period, Id_survey))

# Match periods to surveys ----
## The landscape-metrics tables give each period's ssp at 300 m and the data year its layer applies to ('past' and 'ubc' only)
Lsm_300 <- map(Lc_periods, \(p) read_excel(list.files(Lsm_dir, pattern = paste0("Lsm_df_", p), full.names = TRUE)) |> filter(buffer == Buffer_m) |> select(Id_survey, ssp, any_of("data_year")) |> mutate(Lc_period = p)) |>
  bind_rows() |>
  inner_join(Ssp_period, by = c("Lc_period", "Id_survey"))
Lsm_l <- split(Lsm_300, Lsm_300$Lc_period) |> map(\(d) select(d, -Lc_period) |> rename(any_of(c(Year = "data_year"))))
Cw_cols <- c("ssp", "Ssp_area_poly", "Ssp_cw", "Crown_pct")
Event_covs_pcs <- read_csv(paste0(Wrangling_dir, "Excels/Event_covs_pcs.csv"), show_col_types = FALSE)
Event_ssp <- reduce(Lsm_l[c("ubc", "past")], \(acc, d) rows_update(acc, d, by = c("Id_survey", "Year"), unmatched = "ignore"),
                    .init = Event_covs_pcs |> select(Id_survey, Id_survey_no_dc, Year_grp, Year) |> left_join(select(Lsm_l$middle, -Year), by = "Id_survey")) |>
  group_by(Id_survey) |>
  fill(all_of(Cw_cols)) |>
  ungroup() |>
  group_by(Id_survey_no_dc) |>
  arrange(Year, .by_group = TRUE) |>
  fill(all_of(Cw_cols)) |>
  ungroup() |>
  distinct(Id_survey, Year_grp, across(all_of(Cw_cols))) |>
  filter(!is.na(ssp))

## Check: the area-based ssp carried through these steps equals the one the models use, one value per point and year-group
Lsm_used <- read_csv("../Farm_diversity/Derived/Excels/Event_covs_lsm.csv", show_col_types = FALSE) |> distinct(Id_survey, Year_grp, ssp)
Ssp_check <- inner_join(Event_ssp, Lsm_used, by = c("Id_survey", "Year_grp"), suffix = c("", "_used"), relationship = "one-to-one")
stopifnot(nrow(Ssp_check) == nrow(Event_ssp), isTRUE(all.equal(Ssp_check$ssp, Ssp_check$ssp_used)))
cat("Area-based ssp from the polygons vs the landscape-metrics (raster) ssp: r =", round(cor(Event_ssp$ssp, Event_ssp$Ssp_area_poly), 4), ", max abs difference", round(max(abs(Event_ssp$ssp - Event_ssp$Ssp_area_poly)), 2), "points\n")

# Isolated trees within 300 m ----
## One count per point count, as the tree layer has no date
Tree_counts <- Pc_locs |>
  mutate(Tree_n = lengths(st_is_within_distance(geom, Trees, dist = Buffer_m))) |>
  st_drop_geometry() |>
  as_tibble() |>
  select(Id_survey, Id_group, Institution_name, Tree_n) |>
  arrange(Id_survey)
write_csv(Tree_counts, Trees_csv)

# Export ----
Event_ssp |>
  left_join(select(Tree_counts, Id_survey, Tree_n), by = "Id_survey", relationship = "many-to-one") |>
  select(Id_survey, Year_grp, ssp, Ssp_area_poly, Ssp_cw, Crown_pct, Tree_n) |>
  mutate(across(c(Ssp_area_poly, Ssp_cw, Crown_pct), \(x) round(x, 3))) |>
  arrange(Id_survey, Year_grp) |>
  write_csv(Out_csv)

# Figures ----
## Canopy-weighted vs area-based ssp, one point per point count and year-group
Cw_plot <- ggplot(Event_ssp, aes(ssp, Ssp_cw)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.4) +
  labs(x = "Silvopasture cover within 300 m, area-based (%)", y = "Silvopasture canopy cover within 300 m (%)")
ggsave(file.path(Fig_dir, "Ssp_canopy_vs_area.png"), Cw_plot, width = 6, height = 5.5, bg = "white")
print(Cw_plot)

## Isolated-tree count by farm, one value per location (farms surveyed by several institutions carry one Id_group per institution, so the farm is Id_group without its institution prefix), raw and per 10 ha of pasture within 300 m (pasture share averaged over the location's year-groups), to look for farms digitized differently: trees are left standing in pasture, so a farm with much pasture but few digitized trees suggests trees drawn as stands or not drawn
Buffer_ha <- pi * Buffer_m^2 / 1e4
Tree_locs <- read_csv("../Farm_diversity/Derived/Excels/Event_covs_lsm.csv", show_col_types = FALSE) |>
  distinct(Id_survey, Year_grp, intpast) |>
  summarise(Pasture_pct = mean(intpast), .by = Id_survey) |>
  right_join(Tree_counts, by = "Id_survey") |>
  ## Dispersed-tree stands (crowns excluded) averaged over the land-cover periods: share of the 300 m buffer and number of stands < 500 m2, the other way a digitizer could have drawn isolated trees
  left_join(summarise(Ssp_period, Stand_pct = mean(Stand_pct), N_small_stands = mean(N_small_stands), .by = Id_survey), by = "Id_survey") |>
  left_join(distinct(read_csv(paste0(Wrangling_dir, "Excels/Site_covs.csv"), show_col_types = FALSE), Id_survey_no_dc, Habitat), by = join_by(Id_survey == Id_survey_no_dc)) |>
  mutate(Location = str_remove(Id_survey, "^[A-Z]+-"), Farm = str_remove(Id_group, "^[A-Z]+-")) |>
  filter(!is.na(Pasture_pct)) |>
  summarise(Tree_n = first(Tree_n), Pasture_pct = mean(Pasture_pct), Stand_pct = mean(Stand_pct, na.rm = TRUE), N_small_stands = mean(N_small_stands, na.rm = TRUE), .by = c(Farm, Location)) |>
  mutate(Trees_per_10ha_pasture = if_else(Pasture_pct > 0, Tree_n / (Buffer_ha * Pasture_pct / 100) * 10, NA_real_),
         Program = if_else(str_starts(Farm, "AD"), "Distance farms (GAICA)", "Monitoring farms"))
Tree_farm_plot <- Tree_locs |>
  mutate(Farm = fct_reorder(Farm, Tree_n, median)) |>
  pivot_longer(c(Tree_n, Trees_per_10ha_pasture), names_to = "Measure", values_to = "Value") |>
  mutate(Measure = recode(Measure, Tree_n = "Digitized isolated trees within 300 m", Trees_per_10ha_pasture = "Per 10 ha of pasture within 300 m")) |>
  ggplot(aes(Value, Farm, colour = Program)) +
  geom_boxplot(outlier.shape = NA, colour = "grey60") +
  geom_jitter(height = 0.15, width = 0, alpha = 0.7) +
  facet_wrap(~Measure, scales = "free_x") +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme(legend.position = "top")
ggsave(file.path(Fig_dir, "Isolated_trees_by_farm.png"), Tree_farm_plot, width = 10, height = 7, bg = "white")
print(Tree_farm_plot)

# Report ----
cat("\nIsolated trees within 300 m, by farm (one value per location):\n")
Tree_locs |>
  summarise(Trees_per_10ha_pasture = round(10 * sum(Tree_n) / sum(Buffer_ha * Pasture_pct / 100), 1), N_locs = n(), Median = median(Tree_n), Max = max(Tree_n), Share_zero = round(mean(Tree_n == 0), 2), Pasture_pct = round(median(Pasture_pct)), Stand_pct = round(median(Stand_pct), 1), Small_stands = round(median(N_small_stands), 1), .by = c(Program, Farm)) |>
  arrange(Trees_per_10ha_pasture) |>
  print(n = Inf, width = Inf)
cat("\nCanopy-weighted vs area-based ssp: r =", round(cor(Event_ssp$ssp, Event_ssp$Ssp_cw), 3), "; mean", round(mean(Event_ssp$Ssp_cw), 2), "vs", round(mean(Event_ssp$ssp), 2), "%\n")
