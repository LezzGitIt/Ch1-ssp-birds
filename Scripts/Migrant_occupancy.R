## PhD birds in silvopastoral landscapes ##
# Occupancy of Nearctic-Neotropical migrants in Meta during the wintering season: single-species multi-season occupancy models for a focal migrant (Setophaga striata by default), with the primary model's occupancy structure (hinge distance, mean survey date) (the report's base models replace day of season with an arrival-period indicator, 2026-11-16), and the primary detection structure without Prob_breed (a resident breeding-phenology covariate)
# Model family: Bernoulli occupancy and detection with logit links (stPGOcc: spatial intercept, NNGP exponential, 5 neighbours, effective-range prior set per fit variant (free, small, mid-range or extent-wide), no AR(1)), the single-species counterpart of the primary svcTMsPGOcc fit; seed 20261107
# Visits outside the wintering season (before 15 October or after March; first Meta detection 8 October 2019, none in 439 visits from 23 September to 6 October 2019; no April surveys in Meta) are set to NA in y and the detection covariates, so occupancy means use of the location during the non-breeding season; the report's single-species fits use species-specific windows (swin flag) and its habitat-by-rainfall fits drop pasture (nopast flag), see Configuration
# Collector team enters detection as fixed effects (dummies), not a random effect: only 4 teams remain in season, too few to estimate a variance, and detection is identified from the repeat visits within each location-period even where a team covers a whole period
# Survey_date is a control, not an estimable trend: 2016-17 surveyed almost only forest, and silvopasture and pasture have only November 2019 and January-February 2026, which also differ in time of season
# Also fits a multi-species model of the migrants (svcTMsPGOcc, 1 factor, same data and season window; tMsPGOcc with a location random effect in the non-spatial variants) on the Nearctic-Neotropical migrants with at least 2 in-season detections, taken from the 278-species list (same locations, periods and covariates as the 98-species list)
# Scale of effect of landscape canopy cover: one fit per ring radius (ring cover from Scale_of_effect.R), compared by WAIC against the fit without a ring (ring_compare section)
# Run with: caffeinate -is Rscript Scripts/Migrant_occupancy.R [species=<Genus_species>] [swin | sep24] [nopast] <smoke> <tot_free | tot_small | ... | add_free | dist_free | ms_add_free | fref_r<radius>_free_short | ring_compare> > <log> 2>&1

# Load libraries ----
library(tidyverse)
source("Scripts/spOcc_fns.R")

# Configuration ----
run_sections <- commandArgs(trailingOnly = TRUE)
## Focal species for the single-species fits and plots: species=<Genus_species> on the command line (any species in the 98-species list), S. striata by default
Species <- str_remove(c(grep("^species=", run_sections, value = TRUE), "species=Setophaga_striata")[1], "^species=")
tmsOcc_l_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_98_hab_fix.rds"
## Nearctic-Neotropical migrants for the multi-species model (AVONET long-distance migrants breeding in North America with at least 2 in-season detections in Meta), from the 278-species list
Migrants_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_278_hab_fix.rds"
Migrants <- c("Setophaga_striata", "Setophaga_ruticilla", "Piranga_rubra", "Catharus_ustulatus", "Contopus_virens", "Mniotilta_varia", "Piranga_olivacea", "Parkesia_noveboracensis", "Setophaga_fusca", "Setophaga_cerulea")
Scaling_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm_scaling.csv"
Stacked_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm.rds"
Time_levels <- c("13-14.Early", "13-14.Late", "16-17.Early", "16-17.Late", "19.Early", "22.Early", "24.Early", "25-26.Early")
Out_dir <- "Derived/Rdata/Model_output/Temporal/Single_species"
## Scale of effect of landscape canopy cover (user, 2026-11-15): the ring radii, ring table and raster years of the 98-species scale-of-effect fits (Scale_of_effect.R), one forest-reference fit per radius per species, compared by WAIC against the fit without a ring (fref_free_short)
Ring_radii_m <- c(150, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000, 20000)
Rings_csv <- "Derived/Excels/Canopy_rings_by_location.csv"
## 20 km ring (user, 2026-11-15: migrants may perceive the landscape at large scales), in its own table from Scale_of_effect.R
Rings_20km_csv <- "Derived/Excels/Canopy_rings_by_location_20km.csv"
Raster_years <- c(2013, 2014, 2016, 2017, 2019, 2022, 2024)
Ring_species <- c("Setophaga_striata", "Setophaga_ruticilla", "Piranga_rubra")

## Wintering season as Julian days (non-leap calendar; leap years shift by one day): 15 October (288) to 31 March (90); day of season counts from 15 October so the season does not wrap at the new year
Season_start <- 288
Season_end <- 90
## Passage variants of the multi-species model (user, 2026-11-16): 'sep24' on the command line starts the season on 24 September (Julian 267), the first migrant detection in the modelled locations and the start of GAICA's 2019 distance-gradient surveys, which hold the early detections of C. virens and C. ustulatus (none of S. striata); fit tags then start with sep24_ instead of oct15_. spOccupancy's multi-species models take the missing-value pattern of y from the first species, so the window cannot differ by species
if ("sep24" %in% run_sections) Season_start <- 267
## End of passage: location-periods with a mean date before 15 October (Julian 288) are flagged as passage (Passage = 1); always 0 with the 15 October season start
Passage_end <- 288
## End of arrival (user, 2026-11-16): 1 December (Julian 335), by when eBird Status and Trends shows S. striata's non-breeding distribution settled in South America (roughly 22 November, certainly by early December); location-periods with a mean date before it are flagged as arrival (Arrival = 1)
Arrival_end <- 335
## End of arrival for the multi-species model (2026-11-16): 1 November (Julian 305), when per-visit detection of any migrant in forest jumped (0.19 in 15-31 October 2016 vs 0.38-0.53 from November to mid-December), consistent with the earlier arrival of most migrants than S. striata; Arrival_nov1 = 1 before it
Arrival_end_comm <- 305
## Alternative end of arrival for the multi-species model (user, 2026-11-16): 15 November (Julian 319); unlike 1 November it puts the early-November 2019 location-periods (all habitats) in the arrival period, so the indicator varies outside forest; Arrival_nov15 = 1 before it
Arrival_end_nov15 <- 319
## Species-specific season windows for the single-species fits (user, 2026-11-18): 'swin' on the command line takes the focal species' start, end of arrival and season end from this table, and fit tags start with swin_ instead of oct15_. Start from our Meta data: 24 September for S. ruticilla and P. rubra (first detected 28 September and 6 October 2019, the same start as the multi-species model), 15 October for S. striata (as before; first detected at modelled locations 10 November 2019). End of arrival and season end from eBird Status and Trends (ebirdst_runs, 2023 estimates): the species' non-breeding start and non-breeding end; S. ruticilla and P. rubra have no Meta detections after 23 and 21 February (90 March visits), S. striata's last is 29 March
Species_windows <- tribble(
  ~Species, ~Start, ~Arrival_end, ~End,
  "Setophaga_striata", 288, 347, 95,
  "Setophaga_ruticilla", 267, 333, 74,
  "Piranga_rubra", 267, 319, 60
)
if ("swin" %in% run_sections) {
  Win <- Species_windows[Species_windows$Species == Species, ]
  stopifnot(nrow(Win) == 1)
  Season_start <- Win$Start; Arrival_end <- Win$Arrival_end; Season_end <- Win$End
}
## Pasture dropped from the habitat-by-rainfall fits (user, 2026-11-18: a pasture rainfall slope would rest on a handful of detections; pasture stays in the base models, where it carries the silvopasture vs pasture contrast and part of the shared slopes): 'nopast' on the command line sets pasture visits to NA, so pasture locations drop out and habitat is forest (reference) vs silvopasture; fit tags gain a nopast_ prefix
Season_tag <- paste0(if ("nopast" %in% run_sections) "nopast_", if ("swin" %in% run_sections) "swin" else if ("sep24" %in% run_sections) "sep24" else "oct15")

## Day of season enters occupancy as the location-period mean (occupancy is constant within a location-period, whose visits are at most ~2 weeks apart), as a quadratic so occupancy can peak mid-season
## Rainfall: annual (Tot_prec, as the multi-species models) or driest-month rainfall (Driest_prec, WorldClim monthly normals; Morales-Rozo et al. 2025 found it supported, while December-March rainfall correlates 0.96 with annual rainfall in Meta)
Occ_formulas <- list(
  tot = ~ Habitat * Tot_prec + Dist_open + Dist_inside + ssp + Survey_date + Day_season + I(Day_season^2),
  dry = ~ Habitat * Driest_prec + Dist_open + Dist_inside + ssp + Survey_date + Day_season + I(Day_season^2),
  # Additive habitat and rainfall (user, 2026-11-14: neither rainfall interaction is supported, as rainfall varies only among 4-5 farms per habitat); dist adds silvopasture's difference from the pasture slope of distance from the forest edge (Ssp_Dist_open, built by hand as in spOccupancy_temporal.R), which varies within farms
  add = ~ Habitat + Tot_prec + Dist_open + Dist_inside + ssp + Survey_date + Day_season + I(Day_season^2),
  dist = ~ Habitat + Tot_prec + Dist_open + Ssp_Dist_open + Dist_inside + ssp + Survey_date + Day_season + I(Day_season^2),
  # Time-control sensitivity (2026-11-14): Survey_date is confounded with Day_season outside forest (r = 0.96, as pasture and silvopasture were surveyed only in November 2019 and January-February 2026); nodate drops it, period replaces it with survey-period intercepts (2025-26 reference) so habitat contrasts come only from comparisons within periods, and period_shared is period on the two periods surveyed in all three habitats
  nodate = ~ Habitat + Tot_prec + Dist_open + Dist_inside + ssp + Day_season + I(Day_season^2),
  period = ~ Habitat + Tot_prec + Dist_open + Dist_inside + ssp + Period_16_17_Early + Period_16_17_Late + Period_19_Early + Day_season + I(Day_season^2),
  period_shared = ~ Habitat + Tot_prec + Dist_open + Dist_inside + ssp + Period_19_Early + Day_season + I(Day_season^2),
  # Forest-reference variants (2026-11-15): forest is the reference habitat (Habitat_fref), because with pasture as reference the multi-species model pulled each species' pasture level toward a community intercept set by six rare forest-only migrants, inflating every contrast against pasture; the two 2016-17 periods share one intercept (Period_16_17), since 2016-17 late surveyed only forest and its intercept was not identified for species never detected there; day of season is linear (within-period spread is too small to identify a curve once periods have intercepts)
  fref = ~ Habitat_fref + Tot_prec + Dist_open + Dist_inside + ssp + Period_16_17 + Period_19_Early + Day_season,
  # fref_dist adds silvopasture's difference from the pasture slope of distance outside forest; fref_split splits silvopasture into live fence and dispersed trees
  fref_dist = ~ Habitat_fref + Tot_prec + Dist_open + Ssp_Dist_open + Dist_inside + ssp + Period_16_17 + Period_19_Early + Day_season,
  fref_split = ~ Habitat_split + Tot_prec + Dist_open + Dist_inside + ssp + Period_16_17 + Period_19_Early + Day_season,
  # Coding check: fref with pasture as the reference, to separate the effect of the reference level (through the coefficient priors) from that of the merged periods and linear season
  pref = ~ Habitat + Tot_prec + Dist_open + Dist_inside + ssp + Period_16_17 + Period_19_Early + Day_season,
  # fref plus forest cover within 300 m (PLAND forest, 2026-11-15), so the silvopasture-cover and habitat effects are net of nearby forest
  fref_forest = ~ Habitat_fref + Tot_prec + Dist_open + Dist_inside + ssp + Forest_300 + Period_16_17 + Period_19_Early + Day_season,
  # fref with a habitat-by-annual-rainfall interaction (2026-11-15), as in the 98-species primary model
  fref_rain = ~ Habitat_fref * Tot_prec + Dist_open + Dist_inside + ssp + Period_16_17 + Period_19_Early + Day_season
)
## fref plus landscape canopy cover in one ring (fref_r<radius>), every other term unchanged (user, 2026-11-15, as in the 98-species fits)
Occ_formulas <- c(Occ_formulas, map(set_names(Ring_radii_m, paste0("fref_r", Ring_radii_m)), \(radius) update(Occ_formulas$fref, as.formula(paste0("~ . + Canopy_r", radius)))))
## The base models without day of season (user, 2026-11-16, doubts the day-of-season curve): with period intercepts, day of season varies within periods mainly between October-December 2016 and March 2017 (forest only), which share the 2016-17 intercept
Occ_formulas <- c(Occ_formulas, map(c(fref_noday = "fref", fref_r3000_noday = "fref_r3000", fref_r5000_noday = "fref_r5000"), \(occ) update(Occ_formulas[[occ]], ~ . - Day_season)))
## The same with an arrival-period indicator in place of day of season (2026-11-16): within periods it contrasts October-November 2016 with early December 2016 and March 2017 (forest only), as all November 2019 location-periods fall before 1 December
Occ_formulas <- c(Occ_formulas, map(c(fref_arr = "fref_noday", fref_r3000_arr = "fref_r3000_noday", fref_r5000_arr = "fref_r5000_noday"), \(occ) update(Occ_formulas[[occ]], ~ . + Arrival)))
Occ_formulas$fref_r3000_arrnov <- update(Occ_formulas$fref_r3000_noday, ~ . + Arrival_nov1)
Occ_formulas$fref_r3000_arrnov15 <- update(Occ_formulas$fref_r3000_noday, ~ . + Arrival_nov15)
## Habitat-by-rainfall interaction on each current base model (2026-11-18, run with nopast: a pasture rainfall slope would rest on a handful of detections): S. ruticilla and P. rubra, S. striata with its 5 km canopy ring, and the multi-species model
Occ_formulas$fref_rain_arr <- update(Occ_formulas$fref_rain, ~ . - Day_season + Arrival)
Occ_formulas$fref_r5000_rain_arr <- update(Occ_formulas$fref_r5000_arr, ~ . + Habitat_fref:Tot_prec)
## Passage indicator on top of the arrival indicator (nested coding: Arrival covers passage too, so the Passage coefficient is passage relative to arrival), for the sep24 season; arrival ending 1 December or 15 November
Occ_formulas <- c(Occ_formulas, map(c(fref_r3000_arr_pass = "fref_r3000_arr", fref_r3000_arrnov15_pass = "fref_r3000_arrnov15"), \(occ) update(Occ_formulas[[occ]], ~ . + Passage)))
Occ_formulas$fref_r3000_rain_arrnov15_pass <- update(Occ_formulas$fref_r3000_arrnov15_pass, ~ . + Habitat_fref:Tot_prec)
## Habitat-by-cover interactions on the community base model, one per fit (user, 2026-11-19, run with sep24 nopast, so silvopasture's slope differs from forest's): silvopasture cover within 300 m, and canopy cover in the community's 3 km ring; without pasture the distance hinge already gives silvopasture and forest their own distance slopes
Occ_formulas$fref_r3000_ssp_arrnov15_pass <- update(Occ_formulas$fref_r3000_arrnov15_pass, ~ . + Habitat_fref:ssp)
Occ_formulas$fref_r3000_canopy_arrnov15_pass <- update(Occ_formulas$fref_r3000_arrnov15_pass, ~ . + Habitat_fref:Canopy_r3000)
## Detection: team fixed effects, and (team_day) a team-specific Sampling_day_log slope, since some teams may need no warm-up period
Det_formulas <- list(
  team = ~ Forest_bin + Pc_start + I(Pc_start^2) + Sampling_day_log + Team_4 + Team_7 + Team_8,
  team_day = ~ Forest_bin + Pc_start + I(Pc_start^2) + Sampling_day_log * (Team_4 + Team_7 + Team_8),
  ## The sep24 season adds GAICA's 2019 distance-gradient teams (1 and 2), which the 15 October season does not reach
  team_day_sep24 = ~ Forest_bin + Pc_start + I(Pc_start^2) + Sampling_day_log * (Team_1 + Team_2 + Team_4 + Team_7 + Team_8)
)
## Periods surveyed in all three habitats, for the shared-periods variant (drops the forest-only 2016-17 surveys)
Shared_periods <- c("19.Early", "25-26.Early")
## Reference team (the one with most in-season visits, 2025-26) and the visits dropped because their observer was not recorded (the 3 Unillanos visits of November 2019 in the pooled 'Unk' team, no detections)
Team_ref <- "5"
Team_drop <- "10"
Prec_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Prec_df.csv"
## Silvopasture sub-type of each location (live fence or dispersed trees), from the site table, matched through the landscape-metrics table's survey IDs
Site_covs_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Site_covs.csv"
Lsm_path <- "../Farm_diversity/Derived/Excels/Event_covs_lsm.csv"
Seed <- 20261107
Batch_length <- 25
Mcmc <- list(n.batch = 12000, n.burn = 150000, n.thin = 150, n.chains = 3)
## Shorter chains (user, 2026-11-15: convergence is settled for this structure): 120,000 iterations with 60,000 burn-in, thinned by 60, so each chain still keeps 1,000 samples; chosen by a _short suffix on the variant name and tag
Mcmc_short <- list(n.batch = 4800, n.burn = 60000, n.thin = 60, n.chains = 3)
## Intermediate chains for multi-species fits (2026-11-15: 120,000 iterations left the community variances at Rhat 1.12-1.18): 200,000 iterations with 100,000 burn-in, thinned by 100 (1,000 kept per chain); _mid suffix
Mcmc_mid <- list(n.batch = 8000, n.burn = 100000, n.thin = 100, n.chains = 3)
## Doubled chains for multi-species fits (2026-11-16: the among-species variance of the silvopasture coefficient mixed poorly, Rhat 1.12-1.31 across variants): 600,000 iterations with 300,000 burn-in, thinned by 300 (1,000 kept per chain); _long suffix
Mcmc_long <- list(n.batch = 24000, n.burn = 300000, n.thin = 300, n.chains = 3)
## Prior sensitivity (2026-11-15): variance of the normal prior on the occupancy coefficients (single-species) or on their community-level means (multi-species); spOccupancy's default is 2.72; _bvar<value> suffix

## Spatial range priors (phi = 3 / effective range in m): free from 0.5 km (within-farm point spacing) to 66 km (the largest distance between points; a longer range would act as an intercept); small 0.5-3 km; mid 25-30 km; large 60-66 km (extent-wide)
## Among-species variances (multi-species only): inverse-gamma(1, 1) on every occupancy and detection coefficient (2026-11-14), weakly informative (among-species SD median 1.2, 97.5% quantile 6.3, logit scale); the default IG(0.1, 0.1) let the silvopasture variance wander (median 30, Rhat 1.18), as 6 of the 10 migrants were never detected in silvopasture
Tau_ig_migrants <- list(a = 1, b = 1)
Range_priors <- list(free = c(0.5e3, 66e3), small = c(0.5e3, 3e3), mid = c(25e3, 30e3), large = c(60e3, 66e3))
## Fit variants, named <occupancy formula>_<range> (e.g. tot_free), chosen on the command line; tot_ / dry_ use the team detection formula, add_ / dist_ the team_day one, and ms_ fits the multi-species model; fits up to 2026-11-14 afternoon (tags prec_hinge_date_dayseason*) started the season on 1 October and had a collector-team random effect
## spatial = FALSE (multi-species only) fits tMsPGOcc with an unstructured location random effect in place of the spatial factor, as the 98-species non-spatial sensitivity (Multi_spp_100 row 85); the range prior is then unused and the tag says nonsp
make_variant <- \(occ, range, det, multi = FALSE, periods = NULL, tau_ig = NULL, short = FALSE, mid = FALSE, long = FALSE, beta_var = NULL, spatial = TRUE) list(occ.formula = if (spatial) Occ_formulas[[occ]] else update(Occ_formulas[[occ]], ~ . + (1 | Site_num)), spatial = spatial, det.formula = Det_formulas[[det]], priors = c(list(phi.unif = rev(3 / Range_priors[[range]])), if (!is.null(beta_var) && !multi) list(beta.normal = list(mean = 0, var = beta_var))), tau_ig = tau_ig, beta_var = beta_var, mcmc = if (short) Mcmc_short else if (mid) Mcmc_mid else if (long) Mcmc_long else Mcmc, multi = multi, periods = periods,
                                                                                        tag = paste0(Season_tag, "_", if (det == "team") "teamfe" else det, "_", occ, "_", if (spatial) range else "nonsp", if (!is.null(tau_ig)) paste0("_ig", tau_ig$a, tau_ig$b), if (!is.null(beta_var)) paste0("_bvar", beta_var), if (short) "_short", if (mid) "_mid", if (long) "_long"))
Fit_variants <- c(
  cross2(c("tot", "dry"), names(Range_priors)) %>% set_names(map_chr(., \(x) paste(x, collapse = "_"))) %>% map(\(x) make_variant(x[[1]], x[[2]], "team")),
  list(add_free = make_variant("add", "free", "team_day"), dist_free = make_variant("dist", "free", "team_day"), ms_add_free = make_variant("add", "free", "team_day", multi = TRUE),
       nodate_free = make_variant("nodate", "free", "team_day"), period_free = make_variant("period", "free", "team_day"), period_shared_free = make_variant("period_shared", "free", "team_day", periods = Shared_periods),
       ms_period_free = make_variant("period", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants),
       fref_free = make_variant("fref", "free", "team_day"), ms_fref_free = make_variant("fref", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants),
       ms_fref_dist_free = make_variant("fref_dist", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants), ms_fref_split_free = make_variant("fref_split", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants),
       pref_free = make_variant("pref", "free", "team_day"),
       # Forest-cover comparison on the shorter chains: the forest-reference base and the base plus forest cover within 300 m, single- and multi-species
       fref_free_short = make_variant("fref", "free", "team_day", short = TRUE), fref_forest_free_short = make_variant("fref_forest", "free", "team_day", short = TRUE),
       ms_fref_free_short = make_variant("fref", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, short = TRUE), ms_fref_forest_free_short = make_variant("fref_forest", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, short = TRUE),
       ms_fref_rain_free_mid = make_variant("fref_rain", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, mid = TRUE),
       fref_rain_free = make_variant("fref_rain", "free", "team_day"),
       fref_free_bvar1 = make_variant("fref", "free", "team_day", beta_var = 1), fref_free_bvar10 = make_variant("fref", "free", "team_day", beta_var = 10),
       ms_fref_free_bvar1 = make_variant("fref", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, beta_var = 1), ms_fref_free_bvar10 = make_variant("fref", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, beta_var = 10)),
  ## Scale of effect on the shorter chains: fref_r<radius>_free_short, baseline fref_free_short
  map(set_names(paste0("fref_r", Ring_radii_m)), \(occ) make_variant(occ, "free", "team_day", short = TRUE)) %>% set_names(\(occ) paste0(occ, "_free_short")),
  ## The same for the multi-species migrant model (ms_fref_r<radius>_free_short, baseline ms_fref_free_short), on the shorter chains of its existing baseline so WAIC compares
  map(set_names(paste0("fref_r", Ring_radii_m)), \(occ) make_variant(occ, "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, short = TRUE)) %>% set_names(\(occ) paste0("ms_", occ, "_free_short")),
  ## Full-length chains for the best community radii (2026-11-15), against the full-length baseline ms_fref_free
  map(set_names(paste0("fref_r", c(3000, 5000))), \(occ) make_variant(occ, "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants)) %>% set_names(\(occ) paste0("ms_", occ, "_free")),
  ## Full-length single-species chains at S. striata's best radius (2026-11-16), the base model for S. striata in the report, against the full-length baseline fref_free
  list(fref_r5000_free = make_variant("fref_r5000", "free", "team_day")),
  ## Base models without day of season (2026-11-16): S. ruticilla and P. rubra (fref_noday), S. striata (fref_r5000_noday), the multi-species model (fref_r3000_noday), full-length chains
  ## Arrival indicator in place of day of season (2026-11-16), the same species and radii
  ## Multi-species model with doubled chains: arrival ending 1 December (as the single-species fits), 1 November or 15 November
  list(ms_fref_r3000_arr_free_long = make_variant("fref_r3000_arr", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE), ms_fref_r3000_arrnov_free_long = make_variant("fref_r3000_arrnov", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE), ms_fref_r3000_arrnov15_free_long = make_variant("fref_r3000_arrnov15", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE)),
  list(fref_arr_free = make_variant("fref_arr", "free", "team_day"), fref_r5000_arr_free = make_variant("fref_r5000_arr", "free", "team_day"), ms_fref_r3000_arr_free = make_variant("fref_r3000_arr", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants)),
  list(fref_noday_free = make_variant("fref_noday", "free", "team_day"), fref_r5000_noday_free = make_variant("fref_r5000_noday", "free", "team_day"), ms_fref_r3000_noday_free = make_variant("fref_r3000_noday", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants)),
  ## Passage tests on the shorter chains (2026-11-16, run with sep24): passage + arrival ending 1 December, passage + arrival ending 15 November, and arrival ending 1 December alone
  list(ms_fref_r3000_arr_pass_free_short = make_variant("fref_r3000_arr_pass", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, short = TRUE), ms_fref_r3000_arrnov15_pass_free_short = make_variant("fref_r3000_arrnov15_pass", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, short = TRUE), ms_fref_r3000_arr_free_short = make_variant("fref_r3000_arr", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, short = TRUE)),
  ## The passage + 15 November arrival test (the only passage coding whose occupancy terms converged, lowest WAIC of the three) with doubled chains, as a candidate community base model
  list(ms_fref_r3000_arrnov15_pass_free_long = make_variant("fref_r3000_arrnov15_pass", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE)),
  ## Single-species fits for the species-specific windows (2026-11-18, run with swin; the habitat-by-rainfall fits also with nopast): S. ruticilla and P. rubra start on 24 September, so their detection formula carries teams 1 and 2 (the _t12 variants)
  list(fref_arr_t12_free = make_variant("fref_arr", "free", "team_day_sep24"), fref_r5000_rain_arr_free = make_variant("fref_r5000_rain_arr", "free", "team_day"), fref_rain_arr_t12_free = make_variant("fref_rain_arr", "free", "team_day_sep24")),
  ## Multi-species habitat-by-rainfall fit on the community base model, doubled chains (2026-11-18, run with sep24 nopast)
  list(ms_fref_r3000_rain_arrnov15_pass_free_long = make_variant("fref_r3000_rain_arrnov15_pass", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE)),
  ## Multi-species habitat-by-cover fits on the community base model, doubled chains (2026-11-19, run with sep24 nopast)
  list(ms_fref_r3000_ssp_arrnov15_pass_free_long = make_variant("fref_r3000_ssp_arrnov15_pass", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE), ms_fref_r3000_canopy_arrnov15_pass_free_long = make_variant("fref_r3000_canopy_arrnov15_pass", "free", "team_day_sep24", multi = TRUE, tau_ig = Tau_ig_migrants, long = TRUE)),
  ## Community at the chosen 3 km radius without the spatial random effect, and its no-ring baseline (2026-11-16, user: location random effect instead), full-length chains
  list(ms_fref_nonsp = make_variant("fref", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, spatial = FALSE), ms_fref_r3000_nonsp = make_variant("fref_r3000", "free", "team_day", multi = TRUE, tau_ig = Tau_ig_migrants, spatial = FALSE))
)

# Build the covariates ----
tmsOcc_l <- readRDS(tmsOcc_l_path)
Scaling <- read_csv(Scaling_path, show_col_types = FALSE)
oc <- tmsOcc_l$occ.covs

## Hinge distance to the forest edge, as in spOccupancy_temporal.R
Dist_sd <- filter(Scaling, Process == "occupancy", Covariate == "Dist_forest")
Dist_m <- oc$Dist_forest * Dist_sd$SD + Dist_sd$Mean
Is_forest <- matrix(oc$Habitat == "Bosque", nrow(Dist_m), ncol(Dist_m))
oc[c("Dist_open", "Dist_inside")] <- hinge_distance(Dist_m, Is_forest, Dist_sd$SD)
oc$Ssp_Dist_open <- ifelse(matrix(oc$Habitat == "Ssp", nrow(Dist_m), ncol(Dist_m)), oc$Dist_open, 0)

## Mean survey date (decimal years) of each location-period, as in spOccupancy_temporal.R
Stacked_oc <- readRDS(Stacked_path)$occ.covs
Site_levels <- sort(unique(Stacked_oc$Id_survey_no_dc))
stopifnot(identical(as.character(Site_levels), dimnames(tmsOcc_l$y)[[2]]))
Stacked_period <- paste(Stacked_oc$Year_grp, Stacked_oc$Season, sep = ".")
Date_raw <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Year_grp = Stacked_oc$Year_grp, Season = Stacked_oc$Season) %>%
  left_join(survey_dates(read_csv("../Ssp-bird-data-wrangling/Derived/Excels/Event_covs.csv", show_col_types = FALSE)) %>% mutate(Id_survey = as.character(Id_survey)),
            by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
  pull(Survey_date)
stopifnot(!anyNA(Date_raw))
oc$Survey_date <- site_period_matrix((Date_raw - mean(Date_raw)) / sd(Date_raw), Stacked_oc$Id_survey_no_dc, Stacked_period, Site_levels, Time_levels, fill = "period_mean")

## Forest cover within 300 m (PLAND forest, %, from the land cover matched to the survey year), laid out site x period as in spOccupancy_temporal.R (unsurveyed cells from the nearest surveyed period) and scaled over all surveyed location-periods
Forest_raw <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Year_grp = Stacked_oc$Year_grp) %>%
  left_join(distinct(read_csv(Lsm_path, show_col_types = FALSE), Id_survey, Year_grp, forest), by = c("Id_survey", "Year_grp"), relationship = "many-to-one") %>%
  pull(forest)
stopifnot(!anyNA(Forest_raw))
Forest_scale <- c(Mean = mean(Forest_raw), SD = sd(Forest_raw))
oc$Forest_300 <- site_period_matrix((Forest_raw - Forest_scale[["Mean"]]) / Forest_scale[["SD"]], Stacked_oc$Id_survey_no_dc, Stacked_period, Site_levels, Time_levels, fill = "nearest")

## Landscape canopy cover in rings (inner 50 m excluded), matched to the survey year as in the 98-species scale-of-effect fits (latest year within a location x year-group x season, snapped to the nearest raster year), standardised over the 459 surveyed location-periods so the coefficients compare with those fits, unsurveyed cells from the nearest surveyed period; one site x period matrix per radius (Canopy_r<radius>)
Event_years <- read_csv("../Ssp-bird-data-wrangling/Derived/Excels/Event_covs.csv", show_col_types = FALSE) %>%
  distinct(Id_survey, Id_no_dc = Id_survey_no_dc, Year_grp, Season, Year) %>%
  slice_max(by = c(Id_survey, Year_grp, Season), order_by = Year, with_ties = FALSE)
## The 20 km ring is extracted by Scale_of_effect.R from the national canopy maps
stopifnot(file.exists(Rings_20km_csv))
Ring_rows <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Year_grp = as.character(Stacked_oc$Year_grp), Season = as.character(Stacked_oc$Season), Site = Stacked_oc$Id_survey_no_dc, Period = Stacked_period) %>%
  left_join(Event_years, by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
  mutate(Raster_year = map_dbl(Year, \(year) Raster_years[which.min(abs(Raster_years - year))])) %>%
  left_join(filter(bind_rows(read_csv(Rings_csv, show_col_types = FALSE), read_csv(Rings_20km_csv, show_col_types = FALSE)), Radius_m %in% Ring_radii_m), by = c("Id_no_dc" = "Id_survey_no_dc", "Raster_year"), relationship = "many-to-many")
stopifnot(!anyNA(Ring_rows$Canopy_cover_ring), nrow(Ring_rows) == nrow(Stacked_oc) * length(Ring_radii_m))
## SD of each radius's ring cover (percentage points), to report coefficients per 10 percentage points
Ring_sd <- Ring_rows %>% summarise(Ring_SD = sd(Canopy_cover_ring), .by = Radius_m)
Ring_covs <- Ring_rows %>%
  mutate(Canopy_ring = as.numeric(scale(Canopy_cover_ring)), .by = Radius_m) %>%
  split(.$Radius_m) %>%
  map(\(d) site_period_matrix(d$Canopy_ring, d$Site, d$Period, Site_levels, Time_levels, fill = "nearest")) %>%
  set_names(\(radius) paste0("Canopy_r", radius))
oc <- c(oc, Ring_covs)

## Driest-month rainfall (mm) of each location: minimum of the 12 WorldClim monthly normals, averaged over the location's collector-specific point coordinates; scaled after the sites are subset below
Driest_raw <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Site = Stacked_oc$Id_survey_no_dc) %>%
  distinct() %>%
  left_join(read_csv(Prec_path, show_col_types = FALSE) %>% mutate(Driest = pmin(!!!syms(paste0("prec_", 1:12)))) %>% select(Id_survey, Driest), by = "Id_survey", relationship = "many-to-one") %>%
  summarise(Driest = mean(Driest), .by = Site) %>%
  arrange(match(Site, Site_levels)) %>%
  pull(Driest)
stopifnot(length(Driest_raw) == length(Site_levels), !anyNA(Driest_raw))

## Habitat with forest as the reference level, and with silvopasture split into live fence and dispersed trees (forest reference too)
oc$Habitat_fref <- relevel(factor(oc$Habitat), "Bosque")
Ssp_sub <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Site = Stacked_oc$Id_survey_no_dc) %>%
  distinct() %>%
  left_join(read_csv(Lsm_path, show_col_types = FALSE) %>% distinct(Id_survey, Id_no_dc = Id_survey_no_dc), by = "Id_survey", relationship = "many-to-one") %>%
  left_join(read_csv(Site_covs_path, show_col_types = FALSE) %>% distinct(Id_no_dc = Id_survey_no_dc, Habitat_sub), by = "Id_no_dc", relationship = "many-to-one") %>%
  distinct(Site, Habitat_sub) %>%
  arrange(match(Site, Site_levels))
stopifnot(identical(Ssp_sub$Site, Site_levels))
oc$Habitat_split <- factor(case_when(oc$Habitat != "Ssp" ~ as.character(oc$Habitat), Ssp_sub$Habitat_sub == "Cerca viva" ~ "Live_fence", Ssp_sub$Habitat_sub == "Arboles dispersos" ~ "Dispersed_trees"), levels = c("Bosque", "Pastizales", "Live_fence", "Dispersed_trees"))
stopifnot(!anyNA(oc$Habitat_split))

# Single-species data, wintering season only ----
Migrants_l <- readRDS(Migrants_path)
stopifnot(identical(dimnames(Migrants_l$y)[2:3], dimnames(tmsOcc_l$y)[2:3]), identical(Migrants_l$y[Species, , , ], tmsOcc_l$y[Species, , , ]), isTRUE(all.equal(Migrants_l$coords, tmsOcc_l$coords)))
y_ms <- Migrants_l$y[Migrants, , , ]
y <- tmsOcc_l$y[Species, , , ]
Julian_sd <- filter(Scaling, Process == "detection", Covariate == "Julian_day")
Julian <- round(tmsOcc_l$det.covs$Julian_day * Julian_sd$SD + Julian_sd$Mean)
In_season <- !is.na(Julian) & (Julian >= Season_start | Julian <= Season_end)
cat("Visits:", sum(!is.na(y)), "| out of season (set to NA):", sum(!is.na(y) & !In_season), "| detections lost:", sum(y[!In_season], na.rm = TRUE), "\n")
Team <- tmsOcc_l$det.covs$Collector_team_num
cat("Visits of the unrecorded-observer team dropped:", sum(In_season & Team %in% Team_drop, na.rm = TRUE), "| detections lost:", sum(y[In_season & Team %in% Team_drop], na.rm = TRUE), "\n")
In_season <- In_season & !(Team %in% Team_drop)
if ("nopast" %in% run_sections) {
  Is_pasture <- array(oc$Habitat == "Pastizales", dim(Julian))
  cat("Pasture visits dropped:", sum(In_season & Is_pasture, na.rm = TRUE), "| detections lost:", sum(y[In_season & Is_pasture], na.rm = TRUE), "\n")
  In_season <- In_season & !Is_pasture
}
y[!In_season] <- NA
y_ms[rep(!In_season, each = length(Migrants))] <- NA
det_covs <- map(tmsOcc_l$det.covs[setdiff(names(tmsOcc_l$det.covs), c("Prob_breed", "Julian_day", "Collector_team_num"))], \(x) replace(x, !In_season, NA))
## One 0/1 detection covariate per team other than the reference
Teams <- sort(unique(Team[In_season]))
cat("Teams in season:", paste(Teams, collapse = ", "), "| reference:", Team_ref, "\n")
det_covs <- c(det_covs, set_names(map(setdiff(Teams, Team_ref), \(k) ifelse(In_season, as.numeric(Team == k), NA)), paste0("Team_", setdiff(Teams, Team_ref))))

## Day of season: location-period mean over in-season visits, scaled over the surveyed location-periods; unsurveyed cells get the period mean (or 0 when a period has no in-season visits), which only affects predictions there
Day_raw <- ifelse(Julian >= Season_start, Julian - Season_start, Julian + 365 - Season_start)
Day_raw[!In_season] <- NA
Day_lp <- apply(Day_raw, 1:2, \(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
Day_mean <- mean(Day_lp, na.rm = TRUE); Day_sd <- sd(Day_lp, na.rm = TRUE)
Day_scaled <- (Day_lp - Day_mean) / Day_sd
Period_mean <- colMeans(Day_scaled, na.rm = TRUE)
oc$Day_season <- ifelse(is.na(Day_scaled), matrix(replace(Period_mean, is.nan(Period_mean), 0), nrow(y), ncol(y), byrow = TRUE), Day_scaled)
cat("Day of season (days since the season start) over surveyed location-periods: mean", round(Day_mean, 1), "SD", round(Day_sd, 1), "range", paste(round(range(Day_lp, na.rm = TRUE)), collapse = "-"), "\n")

## Drop locations and periods left without any in-season visit
Keep_site <- apply(y, 1, \(v) any(!is.na(v)))
Keep_time <- apply(y, 2, \(v) any(!is.na(v)))
cat("Locations kept:", sum(Keep_site), "of", length(Keep_site), "| periods kept:", paste(Time_levels[Keep_time], collapse = ", "), "\n")
subset_sp <- \(x) if (is.null(dim(x))) x[Keep_site] else if (length(dim(x)) == 2) x[Keep_site, Keep_time, drop = FALSE] else x[Keep_site, Keep_time, , drop = FALSE]
ss_l <- list(
  y = y[Keep_site, Keep_time, ],
  occ.covs = c(map(oc[c("Habitat", "Habitat_fref", "Habitat_split", "Tot_prec", "Dist_open", "Ssp_Dist_open", "Dist_inside", "ssp", "Forest_300", "Survey_date", "Day_season", names(Ring_covs))], \(x) if (is.factor(x)) droplevels(subset_sp(x)) else subset_sp(x)), list(Driest_prec = as.numeric(scale(Driest_raw[Keep_site])))),
  det.covs = map(det_covs, subset_sp),
  coords = tmsOcc_l$coords[Keep_site, ]
)
stopifnot(all(map_lgl(ss_l$det.covs, ~ identical(is.na(.x), is.na(ss_l$y)))), !anyNA(unlist(ss_l$occ.covs)))
## Survey-period intercepts (one 0/1 matrix per period other than 2025-26), for the period variants
Periods_kept <- Time_levels[Keep_time]
Period_ref <- "25-26.Early"
ss_l$occ.covs <- c(ss_l$occ.covs, set_names(map(setdiff(Periods_kept, Period_ref), \(k) matrix(as.numeric(Periods_kept == k), sum(Keep_site), length(Periods_kept), byrow = TRUE)), paste0("Period_", str_replace_all(setdiff(Periods_kept, Period_ref), "[-.]", "_"))))
## One intercept for both 2016-17 periods (the forest-reference variants)
ss_l$occ.covs$Period_16_17 <- matrix(as.numeric(str_detect(Periods_kept, "^16-17")), sum(Keep_site), length(Periods_kept), byrow = TRUE)
## Arrival period: location-period mean date before Arrival_end; unsurveyed cells are set to 0 (after arrival), which only affects predictions there
ss_l$occ.covs$Arrival <- replace(1 * (Day_lp[Keep_site, Keep_time] < Arrival_end - Season_start), is.na(Day_lp[Keep_site, Keep_time]), 0)
ss_l$occ.covs$Arrival_nov1 <- replace(1 * (Day_lp[Keep_site, Keep_time] < Arrival_end_comm - Season_start), is.na(Day_lp[Keep_site, Keep_time]), 0)
ss_l$occ.covs$Arrival_nov15 <- replace(1 * (Day_lp[Keep_site, Keep_time] < Arrival_end_nov15 - Season_start), is.na(Day_lp[Keep_site, Keep_time]), 0)
ss_l$occ.covs$Passage <- replace(1 * (Day_lp[Keep_site, Keep_time] < Passage_end - Season_start), is.na(Day_lp[Keep_site, Keep_time]), 0)
Surveyed_kept <- apply(ss_l$y, 1:2, \(v) any(!is.na(v)))
cat("Forest cover within 300 m (%) by habitat, surveyed location-periods:", paste(names(hab_f <- tapply(ss_l$occ.covs$Forest_300[Surveyed_kept] * Forest_scale[["SD"]] + Forest_scale[["Mean"]], matrix(as.character(ss_l$occ.covs$Habitat), nrow(Surveyed_kept), ncol(Surveyed_kept))[Surveyed_kept], mean)), round(hab_f), collapse = ", "), "| r with Dist_inside", round(cor(ss_l$occ.covs$Forest_300[Surveyed_kept], ss_l$occ.covs$Dist_inside[Surveyed_kept]), 2), "| r with ssp", round(cor(ss_l$occ.covs$Forest_300[Surveyed_kept], ss_l$occ.covs$ssp[Surveyed_kept]), 2), "\n")
cat("Silvopasture sub-types (locations kept):", paste(names(table(ss_l$occ.covs$Habitat_split)), table(ss_l$occ.covs$Habitat_split), collapse = ", "), "\n")
## Multi-species list: the same covariates and coordinates with the migrants' detections and their range indicators
ms_l <- c(list(y = y_ms[, Keep_site, Keep_time, , drop = FALSE]), ss_l[c("occ.covs", "det.covs", "coords")], list(range.ind = Migrants_l$range.ind[Migrants, Keep_site, drop = FALSE]))
stopifnot(identical(ms_l$y[Species, , , ], ss_l$y))
cat("Migrants, in-season detections:", paste0(Migrants, " ", apply(ms_l$y, 1, sum, na.rm = TRUE), collapse = "; "), "\n")
cat("Driest-month rainfall (mm): mean", round(mean(Driest_raw[Keep_site]), 1), "SD", round(sd(Driest_raw[Keep_site]), 1), "| r with Tot_prec", round(cor(Driest_raw[Keep_site], as.matrix(ss_l$occ.covs$Tot_prec)[, 1]), 2), "\n")
cat("Correlation of Day_season with Sampling_day_log (visit level):", round(cor(c(Day_raw[Keep_site, Keep_time, ]), c(ss_l$det.covs$Sampling_day_log), use = "complete.obs"), 2), "\n")

## Detections through time: location-periods with a detection out of surveyed location-periods, by period
Det_lp <- apply(ss_l$y, 1:2, \(v) if (all(is.na(v))) NA else any(v == 1, na.rm = TRUE))
tibble(Period = colnames(Det_lp), Surveyed = colSums(!is.na(Det_lp)), Detected = colSums(Det_lp, na.rm = TRUE), Detections = apply(ss_l$y, 2, sum, na.rm = TRUE)) %>%
  mutate(Naive_occ = round(Detected / Surveyed, 2)) %>%
  print()

# Fit ----
## One independent fit per variant named on the command line; smoke runs each named variant (or tot_free) for 10 batches on 1 chain without saving
fit_variant <- function(variant, smoke = FALSE) {
  len <- if (smoke) list(n.batch = 10, n.burn = 0, n.thin = 1, n.chains = 1) else variant$mcmc
  ## Every in-season team other than the reference needs its own detection dummy, or it is silently pooled with the reference team
  stopifnot(all(paste0("Team_", setdiff(Teams, Team_ref)) %in% all.vars(variant$det.formula)))
  set.seed(Seed)
  start <- Sys.time()
  cat("Started", variant$tag, format(start), "\n")
  mod <- if (variant$multi && !variant$spatial) {
    z_inits <- apply(ms_l$y, 1:3, \(v) as.numeric(sum(v, na.rm = TRUE) > 0))
    # Multi-species, non-spatial: tMsPGOcc takes no coordinates or range indicators (no migrant is clipped in Meta), and the location index enters as an occupancy covariate for the random effect
    dat <- list(y = ms_l$y, occ.covs = c(ms_l$occ.covs, list(Site_num = seq_len(dim(ms_l$y)[2]))), det.covs = ms_l$det.covs)
    spOccupancy::tMsPGOcc(occ.formula = variant$occ.formula, det.formula = variant$det.formula, data = dat, inits = list(z = z_inits), priors = c(if (!is.null(variant$tau_ig)) list(tau.sq.beta.ig = variant$tau_ig, tau.sq.alpha.ig = variant$tau_ig), if (!is.null(variant$beta_var)) list(beta.comm.normal = list(mean = 0, var = variant$beta_var))), n.batch = len$n.batch, batch.length = Batch_length, accept.rate = 0.43, ar1 = FALSE, n.omp.threads = 1, verbose = TRUE, n.report = 2000, n.burn = len$n.burn, n.thin = len$n.thin, n.chains = len$n.chains)
  } else if (variant$multi) {
    z_inits <- apply(ms_l$y, 1:3, \(v) as.numeric(sum(v, na.rm = TRUE) > 0))
    # Multi-species: as the primary model's svcTMsPGOcc call, with 1 factor (anchored on S. striata, the first species) for 10 species
    spOccupancy::svcTMsPGOcc(occ.formula = variant$occ.formula, det.formula = variant$det.formula, data = ms_l, inits = list(z = z_inits), priors = c(list(phi.unif = as.list(variant$priors$phi.unif)), if (!is.null(variant$tau_ig)) list(tau.sq.beta.ig = variant$tau_ig, tau.sq.alpha.ig = variant$tau_ig), if (!is.null(variant$beta_var)) list(beta.comm.normal = list(mean = 0, var = variant$beta_var))), tuning = list(phi = 2), svc.cols = 1, cov.model = "exponential", NNGP = TRUE, n.neighbors = 5, n.factors = 1, n.batch = len$n.batch, batch.length = Batch_length, accept.rate = 0.43, ar1 = FALSE, n.omp.threads = 1, verbose = TRUE, n.report = 2000, n.burn = len$n.burn, n.thin = len$n.thin, n.chains = len$n.chains)
  } else {
    # Shared-periods variant: keep only the named periods, then the locations with an in-season visit in them
    dat <- ss_l
    if (!is.null(variant$periods)) {
      keep_t <- Periods_kept %in% variant$periods
      keep_s <- apply(ss_l$y[, keep_t, , drop = FALSE], 1, \(v) any(!is.na(v)))
      sub <- \(x) if (is.null(dim(x))) x[keep_s] else if (length(dim(x)) == 2) x[keep_s, keep_t, drop = FALSE] else x[keep_s, keep_t, , drop = FALSE]
      dat <- list(y = sub(ss_l$y), occ.covs = map(ss_l$occ.covs, sub), det.covs = map(ss_l$det.covs, sub), coords = ss_l$coords[keep_s, ])
      cat("Periods:", paste(Periods_kept[keep_t], collapse = ", "), "| locations:", sum(keep_s), "| visits:", sum(!is.na(dat$y)), "| detections:", sum(dat$y, na.rm = TRUE), "\n")
    }
    z_inits <- apply(dat$y, 1:2, \(v) as.numeric(sum(v, na.rm = TRUE) > 0))
    spOccupancy::stPGOcc(occ.formula = variant$occ.formula, det.formula = variant$det.formula, data = dat, inits = list(z = z_inits), priors = variant$priors, tuning = list(phi = 2), cov.model = "exponential", NNGP = TRUE, n.neighbors = 5, n.batch = len$n.batch, batch.length = Batch_length, accept.rate = 0.43, ar1 = FALSE, n.omp.threads = 1, verbose = TRUE, n.report = 2000, n.burn = len$n.burn, n.thin = len$n.thin, n.chains = len$n.chains)
  }
  cat("Wall time:", format(Sys.time() - start), "\n")
  if (smoke) return(invisible(mod))
  summary(mod)
  print(spOccupancy::waicOcc(mod))
  if (variant$spatial) cat("Effective range (km), median and 95% CrI:", round(quantile(3 / mod$theta.samples[, grep("phi", colnames(mod$theta.samples))[1]], c(.5, .025, .975)) / 1000, 1), "\n")
  dir.create(Out_dir, recursive = TRUE, showWarnings = FALSE)
  save_path <- file.path(Out_dir, paste0("Occ_mod_", if (variant$multi) paste0(length(Migrants), "_boreal_migrants") else Species, "_", variant$tag, "_", Sys.Date(), ".rds"))
  if (file.exists(save_path)) stop("'", save_path, "' already exists")
  saveRDS(mod, save_path); cat("SAVED", save_path, "\n")
  invisible(mod)
}
Run_variants <- intersect(run_sections, names(Fit_variants))
if ("smoke" %in% run_sections) walk(Fit_variants[if (length(Run_variants)) Run_variants else "tot_free"], fit_variant, smoke = TRUE)
if (!"smoke" %in% run_sections) walk(Fit_variants[Run_variants], fit_variant)

# Rainfall interaction plot ----
## Occupancy against annual and driest-month rainfall by habitat (posterior median and 95% CrI, other covariates at their means, distance terms at the forest edge), drawn over each habitat's observed rainfall range; points are farms, since rainfall varies only between farms: the posterior mean of latent occupancy (z) averaged over the farm's surveyed location-periods, sized by location-periods (detection-corrected, unlike the share of location-periods with a detection, which sits below the curves wherever detection < 1)
if ("plot_rain" %in% run_sections) {
  Rain_fits <- c(Annual = "tot_free", `Driest month` = "dry_free")
  Rain_terms <- c(Annual = "Tot_prec", `Driest month` = "Driest_prec")
  Tot_sd <- filter(Scaling, Process == "occupancy", Covariate == "Tot_prec")
  Rain_scale <- list(Annual = c(Tot_sd$Mean, Tot_sd$SD), `Driest month` = c(mean(Driest_raw[Keep_site]), sd(Driest_raw[Keep_site])))
  Hab_labels <- c(Pastizales = "Pasture", Ssp = "Silvopasture", Bosque = "Forest")
  Hab_colors <- c(Pasture = "#eda100", Silvopasture = "#6250d6", Forest = "#008300")
  Farm_name <- Stacked_oc$Farm_name[match(Site_levels[Keep_site], Stacked_oc$Id_survey_no_dc)]
  Surveyed_lp <- apply(ss_l$y, 1:2, \(v) any(!is.na(v)))
  Rain_mods <- map(Rain_fits, \(fit) readRDS(list.files(Out_dir, paste0(Species, "_oct15_teamfe_", fit, "_.*\\.rds$"), full.names = TRUE)))
  Site_tbl <- tibble(Habitat = Hab_labels[as.character(ss_l$occ.covs$Habitat)], Farm = Farm_name, Annual = as.matrix(ss_l$occ.covs$Tot_prec)[, 1] * Tot_sd$SD + Tot_sd$Mean, `Driest month` = Driest_raw[Keep_site], Lp = rowSums(Surveyed_lp))
  Farm_tbl <- imap(Rain_mods, \(mod, rain) {
    z_mean <- apply(mod$z.samples, 2:3, mean)
    Site_tbl %>%
      mutate(Rainfall = rain, mm = .data[[rain]], Z_sum = rowSums(ifelse(Surveyed_lp, z_mean, 0))) %>%
      summarise(mm = mean(mm), Lp = sum(Lp), Z_sum = sum(Z_sum), .by = c(Rainfall, Habitat, Farm))
  }) %>% list_rbind() %>%
    mutate(Occ = Z_sum / Lp)
  Pred_tbl <- imap(Rain_fits, \(fit, rain) {
    b <- Rain_mods[[rain]]$beta.samples
    term <- Rain_terms[[rain]]
    Site_tbl %>%
      summarise(lo = min(.data[[rain]]), hi = max(.data[[rain]]), .by = Habitat) %>%
      mutate(mm = map2(lo, hi, \(lo, hi) seq(lo, hi, length.out = 50))) %>%
      unnest(mm) %>%
      mutate(Rainfall = rain, x = (mm - Rain_scale[[rain]][1]) / Rain_scale[[rain]][2],
             draws = map2(Habitat, x, \(h, x) {
               code <- names(Hab_labels)[Hab_labels == h]
               eta <- b[, "(Intercept)"] + b[, term] * x
               if (code != "Pastizales") eta <- eta + b[, paste0("Habitat", code)] + b[, paste0("Habitat", code, ":", term)] * x
               plogis(eta)
             }),
             med = map_dbl(draws, median), lcl = map_dbl(draws, quantile, .025), ucl = map_dbl(draws, quantile, .975)) %>%
      select(-draws)
  }) %>% list_rbind()
  Rain_plot <- ggplot(Pred_tbl, aes(mm, med, color = Habitat, fill = Habitat)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(data = Farm_tbl, aes(y = Occ, size = Lp), shape = 21, alpha = 0.7, stroke = 0.4, color = "grey20") +
    facet_wrap(~ Rainfall, scales = "free_x") +
    scale_color_manual(values = Hab_colors, aesthetics = c("color", "fill")) +
    scale_size_area(max_size = 6, name = "Location-periods") +
    labs(x = "Rainfall (mm; WorldClim normals)", y = "Occupancy probability (S. striata)", caption = "Lines: posterior median and 95% CrI from the free-range fits; points: farms, posterior mean occupancy of their surveyed location-periods") +
    theme_bw() + theme(legend.position = "top")
  dir.create("Figures/Striata_ss", recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0("Figures/Striata_ss/", Species, "_rainfall_interaction.png"), Rain_plot, width = 9, height = 5, dpi = 200)
  print(Rain_plot)
}

# Distance interaction plot ----
## Occupancy against distance to the forest edge by habitat from the dist_free fit (posterior median and 95% CrI; other covariates at their means, rainfall at the mean), drawn over each habitat's observed distances (negative = depth inside forest); points are locations binned by 50 m within habitat: posterior mean latent occupancy (z) over the bin's surveyed location-periods, sized by location-periods
if ("plot_dist" %in% run_sections) {
  Hab_labels <- c(Pastizales = "Pasture", Ssp = "Silvopasture", Bosque = "Forest")
  Hab_colors <- c(Pasture = "#eda100", Silvopasture = "#6250d6", Forest = "#008300")
  Dist_mod <- readRDS(list.files(Out_dir, paste0(Species, "_oct15_team_day_dist_free_.*\\.rds$"), full.names = TRUE))
  b <- Dist_mod$beta.samples
  z_mean <- apply(Dist_mod$z.samples, 2:3, mean)
  Surveyed_lp <- apply(ss_l$y, 1:2, \(v) any(!is.na(v)))
  ## Signed distance in m per location (mean over its surveyed periods): outside forest positive, inside forest negative
  Signed_m <- (ss_l$occ.covs$Dist_open - ss_l$occ.covs$Dist_inside) * Dist_sd$SD
  Site_tbl <- tibble(Habitat = Hab_labels[as.character(ss_l$occ.covs$Habitat)], m = map_dbl(seq_len(nrow(Signed_m)), \(i) mean(Signed_m[i, Surveyed_lp[i, ]])),
                     Lp = rowSums(Surveyed_lp), Z_sum = rowSums(ifelse(Surveyed_lp, z_mean, 0)))
  Bin_tbl <- Site_tbl %>%
    mutate(Bin = cut_width(m, 50, boundary = 0)) %>%
    summarise(m = mean(m), Lp = sum(Lp), Z_sum = sum(Z_sum), .by = c(Habitat, Bin)) %>%
    mutate(Occ = Z_sum / Lp)
  Pred_tbl <- Site_tbl %>%
    summarise(lo = min(m), hi = max(m), .by = Habitat) %>%
    mutate(m = map2(lo, hi, \(lo, hi) seq(lo, hi, length.out = 60))) %>%
    unnest(m) %>%
    mutate(draws = map2(Habitat, m, \(h, m) {
      d_open <- pmax(m, 0) / Dist_sd$SD; d_inside <- pmax(-m, 0) / Dist_sd$SD
      eta <- b[, "(Intercept)"] + switch(h,
        Pasture = b[, "Dist_open"] * d_open,
        Silvopasture = b[, "HabitatSsp"] + (b[, "Dist_open"] + b[, "Ssp_Dist_open"]) * d_open,
        Forest = b[, "HabitatBosque"] + b[, "Dist_inside"] * d_inside)
      plogis(eta)
    }), med = map_dbl(draws, median), lcl = map_dbl(draws, quantile, .025), ucl = map_dbl(draws, quantile, .975)) %>%
    select(-draws)
  Dist_plot <- ggplot(Pred_tbl, aes(m, med, color = Habitat, fill = Habitat)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(data = Bin_tbl, aes(y = Occ, size = Lp), shape = 21, alpha = 0.7, stroke = 0.4, color = "grey20") +
    scale_color_manual(values = Hab_colors, aesthetics = c("color", "fill")) +
    scale_size_area(max_size = 6, name = "Location-periods") +
    labs(x = "Distance to forest edge (m; negative = inside forest)", y = paste0("Occupancy probability (", str_replace(Species, "_", " "), ")"), caption = "Lines: posterior median and 95% CrI (dist_free fit); points: locations binned by 50 m, posterior mean occupancy of their surveyed location-periods") +
    theme_bw() + theme(legend.position = "top", legend.box = "vertical", plot.caption = element_text(size = 7))
  dir.create("Figures/Striata_ss", recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0("Figures/Striata_ss/", Species, "_distance_interaction.png"), Dist_plot, width = 7.5, height = 5.5, dpi = 200)
  print(Dist_plot)
}

# Scale-of-effect comparison ----
## Per species (single-species fits) and for the migrant community (multi-species fits, Species = "Community"), WAIC of each ring fit minus that of the fit without a ring on the same chain length (short: fref_free_short / ms_fref_free_short; full: ms_fref_free; below 0 means the ring improves on it), the canopy coefficient (community-level mean for the multi-species fits) per 10 percentage points of ring cover (logit scale, posterior median and 95% CrI) and its P(> 0), the largest Rhat, and the spatial random effect: effective range (km), and for the single-species fits its variance and the SD across locations of its posterior mean (how much spatial pattern it carries). Run with ring_compare
### Multi-species fits have one spatial factor whose first species' (S. striata's) loading is fixed at 1 for identifiability, so each species' spatial variance is its loading squared, and S. striata's is fixed at 1; the loadings squared of Ring_species go to their own table
if ("ring_compare" %in% run_sections) {
  find_fit <- function(species, variant) {
    path <- list.files(Out_dir, paste0("^Occ_mod_", species, "_", Fit_variants[[variant]]$tag, "_\\d{4}-\\d{2}-\\d{2}\\.rds$"), full.names = TRUE)
    stopifnot(length(path) <= 1)
    if (length(path) == 1) path else NA_character_
  }
  q3 <- \(x, name) set_names(as.list(quantile(x, c(0.5, 0.025, 0.975))), paste0(name, c("", "_lo", "_hi")))
  ring_stats <- function(path, radius) {
    mod <- readRDS(path)
    range_km <- 3 / mod$theta.samples[, grep("phi", colnames(mod$theta.samples))[1]] / 1000
    out <- tibble(WAIC = spOccupancy::waicOcc(mod)[["WAIC"]], Max_Rhat = max(unlist(mod$rhat), na.rm = TRUE), !!!q3(range_km, "Range_km"))
    out <- if (is.null(mod$lambda.samples)) {
      mutate(out, !!!q3(mod$theta.samples[, "sigma.sq"], "Spatial_var"), SD_w_sites = sd(colMeans(mod$w.samples)))
    } else {
      mutate(out, Loadings_sq = list(map(set_names(Ring_species), \(sp) q3(mod$lambda.samples[, grep(paste0("^", sp, "-"), colnames(mod$lambda.samples))[1]]^2, "Spatial_var")) %>% map(as_tibble) %>% list_rbind(names_to = "Species_loading")))
    }
    if (is.na(radius)) return(out)
    b <- (if (is.null(mod$beta.comm.samples)) mod$beta.samples else mod$beta.comm.samples)[, paste0("Canopy_r", radius)] / Ring_sd$Ring_SD[Ring_sd$Radius_m == radius] * 10
    mutate(out, Canopy_median = median(b), Canopy_lo = quantile(b, 0.025), Canopy_hi = quantile(b, 0.975), P_pos = mean(b > 0))
  }
  Ring_fits <- bind_rows(expand_grid(Species = c(Ring_species, "Community"), Chains = "short", Radius_m = c(NA, Ring_radii_m)),
                         expand_grid(Species = "Community", Chains = "full", Radius_m = c(NA, 3000, 5000))) %>%
    mutate(Variant = paste0(if_else(Species == "Community", "ms_", ""), if_else(is.na(Radius_m), "fref_free", paste0("fref_r", Radius_m, "_free")), if_else(Chains == "short", "_short", "")),
           Path = map2_chr(if_else(Species == "Community", paste0(length(Migrants), "_boreal_migrants"), Species), Variant, find_fit)) %>%
    filter(!is.na(Path)) %>%
    mutate(Stats = map2(Path, Radius_m, ring_stats)) %>%
    unnest(Stats) %>%
    mutate(dWAIC = WAIC - WAIC[is.na(Radius_m)], .by = c(Species, Chains))
  Ring_compare <- Ring_fits %>%
    select(Species, Chains, Radius_m, WAIC, dWAIC, Canopy_median, Canopy_lo, Canopy_hi, P_pos, Max_Rhat, starts_with("Range_km"), starts_with("Spatial_var"), SD_w_sites) %>%
    mutate(across(where(is.double), \(x) round(x, 3))) %>%
    arrange(Species, Chains, Radius_m)
  write_csv(Ring_compare, "Derived/Excels/Migrant_ss_ring_comparison.csv")
  print(as.data.frame(Ring_compare))
  Ring_fits %>%
    filter(Species == "Community") %>%
    select(Chains, Radius_m, Loadings_sq) %>%
    unnest(Loadings_sq) %>%
    mutate(across(where(is.double), \(x) round(x, 3))) %>%
    arrange(Species_loading, Chains, Radius_m) %>%
    write_csv("Derived/Excels/Migrant_ms_ring_spatial_loadings.csv")

  ## ΔWAIC and the canopy coefficient against ring radius, one colour per species
  Ring_plot <- Ring_compare %>%
    filter(!is.na(Radius_m), Chains == "short") %>%
    pivot_longer(c(dWAIC, Canopy_median), names_to = "name", values_to = "value") %>%
    mutate(lo = if_else(name == "dWAIC", NA_real_, Canopy_lo), hi = if_else(name == "dWAIC", NA_real_, Canopy_hi),
           name = recode(name, dWAIC = "WAIC minus WAIC without a ring (lower is better)", Canopy_median = "Canopy coefficient\n(logit, per 10 percentage points)"),
           Species = str_replace(Species, "_", " ")) %>%
    ggplot(aes(x = Radius_m / 1000, y = value, colour = Species)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.08), na.rm = TRUE) +
    geom_line(position = position_dodge(width = 0.08)) +
    scale_x_log10() +
    facet_wrap(~ name, ncol = 1, scales = "free_y") +
    labs(x = "Outer ring radius (km)", y = NULL, colour = NULL) +
    theme_bw() + theme(legend.position = "top", legend.text = element_text(face = "italic"))
  dir.create("Figures/Striata_ss", recursive = TRUE, showWarnings = FALSE)
  ggsave("Figures/Striata_ss/Migrant_ss_ring_profile.png", Ring_plot, width = 7, height = 6.5, dpi = 200)
  print(Ring_plot)
}
