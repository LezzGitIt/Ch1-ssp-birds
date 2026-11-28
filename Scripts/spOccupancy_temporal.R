## PhD birds in silvopastoral landscapes ##
# Multi-season (temporal) spatial multi-species occupancy model (spOccupancy::svcTMsPGOcc with svc.cols = 1, i.e. a spatially varying intercept only), adapted from Doser's space-time vignette (Extra_scripts/spaceTimeModelsHTML.Rmd): each physical location is one site with one spatial effect, and the year-group x season survey periods are a separate time dimension, so stacked copies of a location no longer need jittered coordinates
# Model family: Bernoulli occupancy and Bernoulli detection with logit links (Polya-Gamma sampler); spatial NNGP exponential process with 2 latent factors, their loadings started at an earlier fit's posterior means; no AR(1) temporal effect (periods are unequally spaced 2013-2026), time enters as a fixed effect
# Distance to the forest edge is hinge-coded: Dist_open (distance outside forest, 0 in forest) and Dist_inside (depth inside forest, 0 outside), so each side of the edge has its own slope and the habitat contrasts sit at the edge
# svcTMsPGOcc with svc.cols = 1 is the same model as stMsPGOcc (seeded runs give bit-identical draws, 2026-11-02) but, unlike stMsPGOcc, accepts range.ind (species x location, so clipping is spatial, not seasonal)
# Run headless with: nohup caffeinate -is Rscript -e 'source("Scripts/spOccupancy_temporal.R")' <model(s)> > <log> 2>&1 &  (models: prec_hinge, dist_hinge, year_hinge, hab_ssp, split_hinge, and the habitat-context and habitat-conversion models ctx1500, ctx1500_noht, full1500 and its variants; smoke for a 10-batch, 1-chain dimension check of every model, smoke_<model> for one); or set Rstudio_sections and source the script in RStudio; source() parses the whole file first, so editing the script during a fit cannot change it

# Load libraries ----
library(tidyverse)
library(spOccupancy)
# Pin functions that other attached packages mask (a plain assignment rather than conflicted, which turns every ambiguous call into an error)
select <- dplyr::select
filter <- dplyr::filter
source("Scripts/spOcc_fns.R")

# Configuration ----
## Species set: "97" (species detected on at least 19 visits) or "275" (all species, range-clipped), both on the land cover as corrected on 2026-11-27 (the fits to date used the earlier land cover: 98 and 278 species, 459 location-periods); both share the 457 surveyed location-periods and every covariate, and fit names carry the species count
Species_set <- "97"
## Stacked list the temporal list is built from (one row per location x year-group x season), with the scaling of its covariates
msOcc_l_path <- c(`97` = "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_97_lc2711.rds", `275` = "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_275_lc2711.rds")[[Species_set]]
tmsOcc_l_path <- c(`97` = "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_97_lc2711.rds", `275` = "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_275_lc2711.rds")[[Species_set]]
Scaling_path <- str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv")
## Farm_diversity's landscape metrics (forest = % forest within 300 m, from the land cover matched to the survey year; survey ID to location code), the wrangling repo's site table (silvopasture type per location) and visit table (survey dates)
Lsm_path <- "../Farm_diversity/Derived/Excels/Event_covs_lsm.csv"
Site_covs_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Site_covs.csv"
Event_covs_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Event_covs.csv"
## Landscape canopy cover and woody canopy height in the ring from 50 m out to Ring_radius_m, per surveyed row (written by the scale-of-effect script, which chose 1.5 km), and the canopy-weighted silvopasture cover within 300 m (written by the silvopasture-canopy script)
Fit_covs_path <- "Derived/Excels/Scale_of_effect_fit_covs.csv"
Ring_radius_m <- 1500
Ssp_canopy_path <- "Derived/Excels/Ssp_canopy_by_survey.csv"

## Site = physical location; time = year-group x season in chronological order (replicates stay the within-season visits, so closure is the same as in the stacked model)
Site_col <- "Id_survey_no_dc"
Time_cols <- c("Year_grp", "Season")
Time_levels <- c("13-14.Early", "13-14.Late", "16-17.Early", "16-17.Late", "19.Early", "22.Early", "24.Early", "25-26.Early")
## The multi-season models need occupancy covariates in every site x period cell, including periods a location was not surveyed (no detection data there): covariates that describe the period are filled with the period mean, all other time-varying covariates with the location's value from its nearest surveyed period
Period_covs <- "Survey_date"

## Models to fit (formula names below, or smoke / smoke_<name>): from the command line, or, when the script is sourced in RStudio, from Rstudio_sections (empty: sourcing builds the data only)
Rstudio_sections <- "add_split"
## Minutes to wait after building the data and running any smoke tests before the full fits start: staggers two fits launched together (e.g. two factor counts in two RStudio sessions) so their end-of-fit memory peaks do not coincide; 0 = start at once
Start_delay_min <- 0
run_sections <- if (interactive()) Rstudio_sections else commandArgs(trailingOnly = TRUE)

## Time covariate: Year1 (first year of the year-group) reproduces the fitted model set (Multi_spp_100 rows 107-109, 111, 126); Survey_date (mean survey date of the location-period; user, 2026-11-09) for refits, whose names gain _date
Time_cov <- "Year1"
## Formulas, one per model; the habitat contrasts and their interactions are the inferential targets, time is a control
occ_formulas <- list(
  # Primary model: habitat-specific response to annual precipitation
  prec_hinge = ~ Habitat * Tot_prec + Dist_open + Dist_inside + ssp + Year1,
  # Habitat-specific slope of distance outside forest: Ssp_Dist_open is silvopasture's difference from the pasture slope, built by hand because Habitat * Dist_open would add an all-zero forest column sampled from its prior alone
  dist_hinge = ~ Habitat + Dist_open + Ssp_Dist_open + Dist_inside + Tot_prec + ssp + Year1,
  # Habitat-specific trend over time
  year_hinge = ~ Habitat * Year1 + Tot_prec + Dist_open + Dist_inside + ssp,
  # Habitat-specific slopes of silvopasture cover within 300 m (spillover into pasture, buffering of forest), at fixed forest cover within 300 m, so more silvopasture means less pasture rather than less forest; in silvopasture the cover includes the point's own patch
  hab_ssp = ~ Habitat * ssp + Forest_300 + Tot_prec + Dist_open + Dist_inside + Year1,
  # Primary model with silvopasture split into live fence and dispersed trees (pasture reference)
  split_hinge = ~ Habitat_split * Tot_prec + Dist_open + Dist_inside + ssp + Year1
) %>%
  map(\(f) as.formula(str_replace_all(deparse1(f), "Year1", Time_cov)))
Hinge_models <- names(occ_formulas)

## Habitat-context and habitat-conversion models (moved from the scale-of-effect script, 2026-11-24; their 98-species fits are in Model_output/Temporal/Scale_of_effect/, named soe_<model>): time enters as Year_cont, the mean survey date read as a long-term trend, and landscape canopy cover and woody canopy height as the 1.5 km ring terms
## Context dependence of habitat on landscape forest (user, 2026-11-21): habitat by ring cover at the chosen 1.5 km radius, with the cover-by-height landscape forest adjustment as primary ('ctx1500') and without height as the sensitivity check ('ctx1500_noht'); habitat and rainfall additive. Data check (2026-11-21): 1.5 km cover overlaps across habitats (10th-90th percentile: silvopasture 12-44%, forest 21-44%, pasture 22-51%), but pasture spans only 17-56% on 6 farms, so silvopasture-vs-pasture contrasts below ~17% cover are extrapolated; no point count has low cover with tall canopy (r 0.68 overall, 0.82 within silvopasture)
## Model for the habitat-conversion predictions ('full1500', user 2026-11-21): every habitat interaction whose community-level CrI excluded 0 when fitted alone, together in one model: habitat by 1.5 km ring cover (ctx1500), habitat by silvopasture cover within 300 m (Multi_spp_100 row 111) and habitat by rainfall (row 107: silvopasture 0.29 (-0.01, 0.58), P(>0) 0.96; kept although its CrI touches 0, user 2026-11-21). Dropped as unsupported: silvopasture's distance-to-edge slope (row 108: 0.08 (-0.19, 0.36))
## Forest cover within 300 m added to the habitat-conversion model (user, 2026-11-22), so the habitat-specific ssp slopes are estimated holding forest within 300 m fixed: at forest points, at fixed pasture, each extra percentage point of ssp comes with 0.42 fewer points of forest, so without it the forest-point ssp slope can absorb forest loss, which the conversion scenarios (silvopasture replacing pasture only) would inherit. 'full1500_f300' adds one forest slope for all habitats; 'full1500_hf300' lets it differ by habitat, which adjusts the forest-point ssp slope even if forest matters more at forest points, at the risk of slow convergence (at pasture points forest within 300 m correlates -0.90 with distance to the forest edge)
## Riparian forest indicator added to 'full1500_f300' (user, 2026-11-22): riparian forest points carry most of the high ssp at forest points (mean 14% vs 7.5% at secondary forest points), so the forest-point ssp slope may be picking up forest type; 'full1500_f300_rip' tests whether it shrinks once riparian forest has its own intercept
## Quadratic landscape cover added to 'full1500_f300' (user, 2026-11-23): one squared term of the standardised 1.5 km ring cover, with a community mean and species-level slopes like every occupancy covariate, so each species' response to landscape cover can peak at a different cover; not interacted with habitat or height, since only ~9-13 independent 1.5 km landscapes inform the landscape terms. A community-level squared coefficient below 0 means occupancy peaks at intermediate landscape cover
## Canopy-weighted silvopasture cover within 300 m in place of the area-based ssp (user, 2026-11-23): area-based ssp counts a sparse dispersed-tree stand like a dense one and overstates dispersed-tree canopy about threefold; 'full1500_f300_cw' is 'full1500_f300' with ssp replaced by Ssp_cw (dispersed-tree stands weighted by their imagery tree cover, live fences, intensive silvopasture and isolated-tree crowns at 100%, forage banks excluded), the model the landscape-silvopasture scenario would simulate from. The isolated-tree count is left out until its consistency across farms is judged (user, 2026-11-23)
Context_base_formula <- ~ Habitat + Tot_prec + Dist_open + Dist_inside + ssp + Year_cont
occ_formulas <- c(occ_formulas, list(
  ctx1500 = update(Context_base_formula, ~ . + Habitat * Cover_h1500 + Cover_h1500 * Height_h1500),
  ctx1500_noht = update(Context_base_formula, ~ . + Habitat * Cover_h1500),
  full1500 = ~ Habitat * (Cover_h1500 + ssp + Tot_prec) + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont,
  full1500_f300 = ~ Habitat * (Cover_h1500 + ssp + Tot_prec) + Forest_300 + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont,
  full1500_hf300 = ~ Habitat * (Cover_h1500 + ssp + Tot_prec + Forest_300) + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont,
  full1500_f300_rip = ~ Habitat * (Cover_h1500 + ssp + Tot_prec) + Forest_300 + Riparian + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont,
  full1500_f300_q = ~ Habitat * (Cover_h1500 + ssp + Tot_prec) + Forest_300 + I(Cover_h1500^2) + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont,
  full1500_f300_cw = ~ Habitat * (Cover_h1500 + Ssp_cw + Tot_prec) + Forest_300 + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont
))
## Additive model with silvopasture split by type, the 278-species base model and the one the number of spatial factors is chosen on (user, 2026-11-24): every term additive, and every covariate the causal diagram requires for the direct effect of local habitat and the effect of silvopasture within 300 m (climate, forest within 300 m, landscape forest, time, the distance to the forest edge, and the spatial factors for unmeasured spatial causes); forest cover within 300 m is held so the silvopasture slope compares silvopasture with pasture rather than with forest, and landscape forest enters as 1.5 km canopy cover by woody canopy height, as in the habitat-conversion models
occ_formulas$add_split <- ~ Habitat_split + Tot_prec + Dist_open + Dist_inside + Ssp_cw + Forest_300 + Cover_h1500 * Height_h1500 + Year_cont
## Habitat-specific trend over time on the DAG base (user, 2026-11-27): kept out of the habitat-conversion models because the time interaction does not change their predictions, but a target in its own right
occ_formulas$add_split_year <- ~ Habitat_split * Year_cont + Tot_prec + Dist_open + Dist_inside + Ssp_cw + Forest_300 + Cover_h1500 * Height_h1500
## Habitat-conversion model with silvopasture split by type (user, 2026-11-27): every habitat interaction of 'full1500_f300_cw' except time, with live fence and dispersed trees each given their own slopes; fitted first, with the unsplit 'full1500_f300_cw' as the fallback if it proves too complex to converge
occ_formulas$full1500_f300_cw_split <- ~ Habitat_split * (Cover_h1500 + Ssp_cw + Tot_prec) + Forest_300 + Cover_h1500 * Height_h1500 + Dist_open + Dist_inside + Year_cont
## Fit name: model, then _date for the hinge models refitted with the mean survey date, the factor count, _inits when the loadings start from the earlier fit and _short for the shorter chains (so the existing 98-species fits keep their names, e.g. prec_hinge_f2_inits)
fit_tag <- \(model) paste0(model, if (Time_cov == "Survey_date" && model %in% Hinge_models) "_date", "_f", N_factors, if (Lambda_start == "earlier_fit") "_inits", if (Chain_length == "short") "_short")
det.formula <- ~ Forest_bin + Pc_start + I(Pc_start^2) + Prob_breed + Sampling_day_log + (1 | Collector_team_num)

## MCMC settings: 60,000 iterations, 30,000 burn-in, thinned by 75 (400 kept per chain); the species x site x period posterior arrays (z, psi, likelihood) are 3.6x the stacked size, and at 3,000 kept samples assembling them exceeded R's 18 GB memory limit (2026-10-30)
batch.length <- 25
mcmc_long <- list(n.batch = 2400, n.burn = 30000, n.thin = 75, n.chains = 3)
## Shorter chains for screening runs such as choosing the number of factors: 24,000 iterations, 12,000 burn-in, thinned by 30 (still 400 kept per chain), 40% of the full length; row 107 refitted at this length matched it (all community coefficients within 0.08 posterior SD, max community Rhat 1.04), as in the scale-of-effect fits; "short" adds _short to the fit name
mcmc_short <- list(n.batch = 960, n.burn = 12000, n.thin = 30, n.chains = 3)
Chain_length <- "long"
n.neighbors <- 5
priors <- list(phi.unif = list(3 / 25e3, 3 / 23e3)) # effective spatial range 23-25 km
ar1 <- FALSE
## Spatial factors: 2 for the 98 species (with 3, factor 2 emptied even when started informative; Multi_spp_100 row 110); 2 for the 278 species (Meta_gradient rows 79-82, 2026-11-26): once mirror-image chains are sign-aligned both factors are informative and converged, and the 3- and 4-factor fits carry the same second axis rotated across their extra factors, with no informative third axis
N_factors <- 2
## Anchor species: spOccupancy fixes the first N_factors species' loadings (1 on their own factor, 0 on later factors) to identify the factors, so they should be common, well-detected species with contrasting habitat use; they move to the front in this order, the rest keep the data list's order, and species beyond N_factors have no effect. Tyrannus melancholicus and Manacus manacus lead the data list (the wrangling script's Order_meta), so the 2-factor fits keep their anchors; Sturnella magna is the pasture species (detected at 43% of pasture location-periods, 3% of forest), Momotus momota and Mesembrinibis cayennensis forest species, Thraupis episcopus a silvopasture-leaning species (user, 2026-11-24)
Anchor_species <- c("Tyrannus_melancholicus", "Manacus_manacus", "Sturnella_magna", "Milvago_chimachima", "Crotophaga_ani", "Momotus_momota", "Thraupis_episcopus", "Mesembrinibis_cayennensis")
## Loadings starting values: "earlier_fit" starts factors 1-2 at an earlier fit's posterior-mean loadings, species that fit lacked and factors beyond 2 at 0; "zero" is spOccupancy's default (all free loadings at 0). spOccupancy uses the same loadings start in every chain, so from "zero" each chain's early path decides the sign of a factor whose anchor loads weakly, and chains can settle in mirror-image modes (the 278-species factor screen, where both anchors load weakly); starting every chain at sign-aligned loadings keeps them in one mode
## 98 species: the 2-factor signed-distance fit (row 90), where factor 2 was informative, because 2-factor runs started at 0 often collapsed it (Jeff Doser's remedy for loadings that won't converge). 278 species: the 2-factor add_split short fit with its chains sign-aligned to the majority mode (Meta_gradient row 82), written on the laptop from the loadings extracted on the remote (2026-11-26), so it has no source fit here
Lambda_start <- "earlier_fit"
Lambda_source_path <- c(`97` = "Derived/Rdata/Model_output/Temporal/Occ_mod_98spp_prec_nofarm_long_f2_2026-11-05.rds", `275` = NA)[[Species_set]]
Lambda_inits_path <- c(`97` = "Derived/Rdata/Model_output/Temporal/lambda_inits_prec_nofarm_long_f2.rds", `275` = "Derived/Rdata/Model_output/Temporal/lambda_inits_278spp_add_split_f2_aligned.rds")[[Species_set]]
## Settings passed to fit_temporal() (spOcc_fns.R); every fit is seeded so a rerun reproduces its draws exactly
fit_settings <- list(det.formula = det.formula, priors = priors, n.neighbors = n.neighbors, n.factors = N_factors, batch.length = batch.length, ar1 = ar1, mcmc = if (Chain_length == "short") mcmc_short else mcmc_long, out_dir = "Derived/Rdata/Model_output/Temporal", seed = 20261107)

# Build the temporal list ----
# Reshape the stacked list (rows = site x time) into the multi-season format of the vignette's hbefTrends data: y is species x site x time x replicate, site-constant occupancy covariates are site vectors, time-varying ones are site x time matrices, detection covariates are site x time x replicate arrays, and coords has one row per site; built once, then read from disk
msOcc_l <- readRDS(msOcc_l_path)
if (file.exists(tmsOcc_l_path)) {
  tmsOcc_l <- readRDS(tmsOcc_l_path)
  # Stop if the saved temporal list was built from a different stacked list (e.g. msOcc_l_path was changed but tmsOcc_l_path was not renamed)
  if (!identical(dimnames(tmsOcc_l$y)[[1]], dimnames(msOcc_l$y)[[1]]) || sum(tmsOcc_l$y, na.rm = TRUE) != sum(msOcc_l$y, na.rm = TRUE)) stop("'", tmsOcc_l_path, "' does not match '", msOcc_l_path, "' (different species or detections); give tmsOcc_l_path a new name so the temporal list is rebuilt")
  if (is.null(tmsOcc_l$range.ind) && any(msOcc_l$range.ind == 0)) message("Saved temporal list has no range.ind (built before clipping was added): spatial fits from it are unclipped")
} else {
  oc <- msOcc_l$occ.covs
  ## Site and period of every stacked row
  Site_levels <- sort(unique(oc[[Site_col]]))
  site_i <- match(oc[[Site_col]], Site_levels)
  time_i <- match(do.call(paste, c(oc[Time_cols], sep = ".")), Time_levels)
  stopifnot(!anyNA(time_i), !anyDuplicated(cbind(site_i, time_i)))
  n_site <- length(Site_levels); n_time <- length(Time_levels); n_rep <- dim(msOcc_l$y)[3]
  ## Detections: species x site x period x replicate, NA where a location was not surveyed in a period
  y <- array(NA_real_, c(dim(msOcc_l$y)[1], n_site, n_time, n_rep), dimnames = list(dimnames(msOcc_l$y)[[1]], Site_levels, Time_levels, NULL))
  for (k in seq_len(n_rep)) for (s in seq_len(dim(y)[1])) y[s, , , k][cbind(site_i, time_i)] <- msOcc_l$y[s, , k]
  ## Occupancy covariates: constant within every site -> site vector, otherwise numeric site x time matrix with unsurveyed cells filled (Period_covs with the period mean, others with the location's value from its nearest surveyed period, ties going to the earlier period); text columns and varying factors (row IDs, the time columns) are dropped
  occ_covs <- oc %>% select(-all_of(c(Site_col, Time_cols))) %>% as.list() %>%
    imap(function(v, nm) {
      if (!is.numeric(v) && !is.factor(v)) return(NULL)
      if (all(tapply(seq_along(v), site_i, function(r) n_distinct(v[r]) == 1))) return(v[match(seq_len(n_site), site_i)])
      if (!is.numeric(v)) return(NULL)
      m <- matrix(NA_real_, n_site, n_time, dimnames = list(Site_levels, Time_levels))
      m[cbind(site_i, time_i)] <- v
      # Only the unsurveyed cells take the period mean, so locations surveyed at a different date from the rest of their period keep their own date
      if (nm %in% Period_covs) return(ifelse(is.na(m), matrix(colMeans(m, na.rm = TRUE), n_site, n_time, byrow = TRUE), m))
      fill_nearest_period(m)
    }) %>% compact()
  ## Detection covariates: stacked row x replicate tables -> site x time x replicate arrays, kept only where the survey happened; Y/N text coded 1/0 because arrays must be numeric, and covariates with no values at all (Pc_duration: constant 10 min, so scaling gave NaN) are dropped
  det_covs <- msOcc_l$det.covs %>% keep(~ !all(is.na(as.matrix(.x)))) %>%
    map(function(tbl) {
      m <- as.matrix(tbl); if (is.character(m)) m <- (m == "Y") * 1
      a <- array(NA_real_, c(n_site, n_time, n_rep), dimnames = list(Site_levels, Time_levels, NULL))
      for (k in seq_len(n_rep)) a[cbind(site_i, time_i, k)] <- m[, k]
      a[is.na(y[1, , , ])] <- NA
      a
    })
  ## One coordinate pair per location
  coords <- apply(msOcc_l$coords, 2, function(col) tapply(col, site_i, mean))
  rownames(coords) <- Site_levels
  ## Range clipping: species x location (the stacked list's range.ind is one column per stacked row, and the range maps do not change between periods)
  range_ind <- t(apply(msOcc_l$range.ind, 1, function(v) tapply(v, site_i, max)))
  dimnames(range_ind) <- list(rownames(msOcc_l$range.ind), Site_levels)
  tmsOcc_l <- list(y = y, occ.covs = occ_covs, det.covs = det_covs, coords = coords, range.ind = range_ind)
  ## Checks: same data as the stacked list, and consistent site order and NA pattern
  stopifnot(
    sum(tmsOcc_l$y, na.rm = TRUE) == sum(msOcc_l$y, na.rm = TRUE),
    sum(!is.na(tmsOcc_l$y)) == sum(!is.na(msOcc_l$y)),
    identical(dimnames(tmsOcc_l$y)[[2]], rownames(tmsOcc_l$coords)),
    all(map_lgl(tmsOcc_l$det.covs, ~ identical(is.na(.x), is.na(tmsOcc_l$y[1, , , ])))),
    all(map_lgl(tmsOcc_l$occ.covs, ~ NROW(.x) == n_site)),
    !anyNA(unlist(tmsOcc_l$occ.covs)),
    # Clipping is the same in every stacked row of a location, and no species is detected where it is clipped
    all(t(apply(msOcc_l$range.ind, 1, function(v) tapply(v, site_i, min))) == range_ind),
    sum(y[array(range_ind, dim(y)) == 0], na.rm = TRUE) == 0
  )
  dir.create(dirname(tmsOcc_l_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(tmsOcc_l, tmsOcc_l_path)
}
# Occupancy covariates of the model set ----
### Built on every run rather than stored in the list, which predates them: hinge-coded distance, the silvopasture slope column, forest cover within 300 m, silvopasture split by type, and (for lists built before 2026-11-09) the mean survey date
oc <- tmsOcc_l$occ.covs
Stacked_oc <- msOcc_l$occ.covs
Site_levels <- sort(unique(Stacked_oc[[Site_col]]))
Stacked_period <- do.call(paste, c(Stacked_oc[Time_cols], sep = "."))
to_site_period <- \(v, fill = "nearest") site_period_matrix(v, Stacked_oc[[Site_col]], Stacked_period, Site_levels, Time_levels, fill)
# The layout must reproduce the list's own ssp matrix
stopifnot(identical(dimnames(tmsOcc_l$y)[[2]], as.character(Site_levels)), isTRUE(all.equal(unname(to_site_period(Stacked_oc$ssp)), unname(oc$ssp))))

## Hinge terms from raw metres, assigned by habitat label (Habitat is one value per location, Dist_forest a location x period matrix, so the label is recycled down the period columns)
Dist_scaling <- filter(read_csv(Scaling_path, show_col_types = FALSE), Process == "occupancy", Covariate == "Dist_forest")
Dist_m <- oc$Dist_forest * Dist_scaling$SD + Dist_scaling$Mean
Is_forest <- matrix(oc$Habitat == "Bosque", nrow(Dist_m), ncol(Dist_m))
oc[c("Dist_open", "Dist_inside")] <- hinge_distance(Dist_m, Is_forest, Dist_scaling$SD)
oc$Ssp_Dist_open <- ifelse(matrix(oc$Habitat == "Ssp", nrow(Dist_m), ncol(Dist_m)), oc$Dist_open, 0)
stopifnot(all(oc$Dist_open[Is_forest] == 0), all(oc$Dist_inside[!Is_forest] == 0))

## Forest cover within 300 m, scaled over the surveyed location-periods like the other covariates
Forest_raw <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Year_grp = Stacked_oc$Year_grp) %>%
  left_join(distinct(read_csv(Lsm_path, show_col_types = FALSE), Id_survey, Year_grp, forest), by = c("Id_survey", "Year_grp"), relationship = "many-to-one") %>%
  pull(forest)
stopifnot(!anyNA(Forest_raw))
oc$Forest_300 <- to_site_period((Forest_raw - mean(Forest_raw)) / sd(Forest_raw))

## Silvopasture split into live fence and dispersed trees, one value per location
Split_rows <- split_ssp(Stacked_oc$Habitat, Stacked_oc$Id_survey, Lsm_path, Site_covs_path)
stopifnot(all(tapply(Split_rows, Stacked_oc[[Site_col]], n_distinct) == 1))
oc$Habitat_split <- Split_rows[match(Site_levels, Stacked_oc[[Site_col]])]
stopifnot(identical(oc$Habitat_split %in% c("Live_fence", "Dispersed_trees"), as.character(oc$Habitat) == "Ssp"))

## Mean survey date in decimal years, as 01_Wrangling_ms.R now builds it, scaled over the surveyed location-periods; unsurveyed cells take the period mean
if (is.null(oc$Survey_date)) {
  Date_raw <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Year_grp = Stacked_oc$Year_grp, Season = Stacked_oc$Season) %>%
    left_join(survey_dates(read_csv(Event_covs_path, show_col_types = FALSE)) %>% mutate(Id_survey = as.character(Id_survey)), by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
    pull(Survey_date)
  stopifnot(!anyNA(Date_raw))
  oc$Survey_date <- to_site_period((Date_raw - mean(Date_raw)) / sd(Date_raw), fill = "period_mean")
}
## The same mean survey date under the name the habitat-context models use (the long-term time trend)
oc$Year_cont <- oc$Survey_date

## Landscape canopy cover and woody canopy height in the 1.5 km ring, scaled over the surveyed location-periods; a ring with no canopy has no woody canopy height and takes the mean (0 once scaled), leaving its zero cover to the cover term
Row_keys <- tibble(Id_survey = as.character(Stacked_oc$Id_survey), Year_grp = as.character(Stacked_oc$Year_grp), Season = as.character(Stacked_oc$Season))
Ring_raw <- Row_keys %>%
  left_join(read_csv(Fit_covs_path, col_types = cols(Id_survey = "c", Year_grp = "c", Season = "c")) %>% filter(Radius_m == Ring_radius_m), by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one")
stopifnot(!anyNA(Ring_raw$Cover_pct))
scale_rows <- \(x) replace_na((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE), 0)
oc[[paste0("Cover_h", Ring_radius_m)]] <- to_site_period(scale_rows(Ring_raw$Cover_pct))
oc[[paste0("Height_h", Ring_radius_m)]] <- to_site_period(scale_rows(Ring_raw$Height_m))

## Riparian forest indicator, one value per location
oc$Riparian <- to_site_period(riparian_forest(Stacked_oc$Habitat, Stacked_oc$Id_survey, Lsm_path, Site_covs_path))

## Canopy-weighted silvopasture cover within 300 m, scaled over the surveyed location-periods
Ssp_cw_raw <- Row_keys %>%
  left_join(read_csv(Ssp_canopy_path, col_types = cols(Id_survey = "c", Year_grp = "c")) %>% select(Id_survey, Year_grp, Ssp_cw), by = c("Id_survey", "Year_grp"), relationship = "many-to-one") %>%
  pull(Ssp_cw)
stopifnot(!anyNA(Ssp_cw_raw))
oc$Ssp_cw <- to_site_period(scale_rows(Ssp_cw_raw))
tmsOcc_l$occ.covs <- oc

## Anchor species first, in the order given (y and range clipping are species-first)
stopifnot(all(Anchor_species %in% dimnames(tmsOcc_l$y)[[1]]), !anyDuplicated(Anchor_species), length(Anchor_species) >= N_factors)
Species_order <- c(Anchor_species, setdiff(dimnames(tmsOcc_l$y)[[1]], Anchor_species))
tmsOcc_l$y <- tmsOcc_l$y[Species_order, , , , drop = FALSE]
if (!is.null(tmsOcc_l$range.ind)) tmsOcc_l$range.ind <- tmsOcc_l$range.ind[Species_order, , drop = FALSE]
str(tmsOcc_l, max.level = 2)

# Loadings starting values ----
## Posterior means of the earlier fit's loadings (species x factor), computed once
if (Lambda_start == "earlier_fit" && !file.exists(Lambda_inits_path)) {
  if (is.na(Lambda_source_path)) stop("Loadings starting values not found: ", Lambda_inits_path)
  Lambda_source <- readRDS(Lambda_source_path)
  saveRDS(matrix(colMeans(Lambda_source$lambda.samples), nrow = length(Lambda_source$sp.names), dimnames = list(Lambda_source$sp.names, NULL)), Lambda_inits_path)
  rm(Lambda_source); gc()
}
## Matched to species by name: the earlier fit's species keep its loadings on factors 1-2, other species and later factors start at 0; spOccupancy requires each anchor's own loading at 1 and its loadings on later factors at 0, and the earlier loadings are relative to its two anchors, so those must stay anchors 1-2 (for 3+ factors the third anchor's loadings on factors 1-2 also start at the earlier fit's values)
lambda_inits <- if (Lambda_start == "earlier_fit") {
  Lambda_earlier <- readRDS(Lambda_inits_path)
  Fit_species <- dimnames(tmsOcc_l$y)[[1]]
  # Species of the earlier fit no longer in the data (below the detection cutoff after a data change) are dropped
  if (!all(rownames(Lambda_earlier) %in% Fit_species)) cat("Loadings starting values: dropped", sum(!rownames(Lambda_earlier) %in% Fit_species), "species not in the data:", setdiff(rownames(Lambda_earlier), Fit_species), "\n")
  Lambda_earlier <- Lambda_earlier[rownames(Lambda_earlier) %in% Fit_species, , drop = FALSE]
  stopifnot(identical(head(Fit_species, 2), head(rownames(Lambda_earlier), 2)))
  Lambda <- matrix(0, length(Fit_species), N_factors, dimnames = list(Fit_species, NULL))
  Lambda[rownames(Lambda_earlier), seq_len(min(2, N_factors))] <- Lambda_earlier[, seq_len(min(2, N_factors))]
  Lambda[upper.tri(Lambda)] <- 0
  diag(Lambda) <- 1
  list(lambda = list(Lambda))
}

# Smoke test ----
# 10 batches on 1 chain, nothing saved: catches formula / dimension errors in about a minute per model; not for judging convergence. 'smoke' checks every model, 'smoke_<model>' one
Smoke_models <- if ("smoke" %in% run_sections) names(occ_formulas) else intersect(str_remove(str_subset(run_sections, "^smoke_"), "^smoke_"), names(occ_formulas))
walk(Smoke_models, \(model) {
  fit_temporal(tmsOcc_l, occ_formulas[[model]], fit_tag(model), settings = fit_settings, inits = lambda_inits, smoke = TRUE)
  cat("Smoke fit OK:", model, "\n")
})

## Staggered start
if (Start_delay_min > 0 && length(setdiff(run_sections, c("smoke", str_subset(run_sections, "^smoke_")))) > 0) {
  cat("Waiting", Start_delay_min, "min before fitting, from", format(Sys.time(), "%H:%M"), "\n")
  Sys.sleep(Start_delay_min * 60)
}

# Primary model: habitat x precipitation ----
if ("prec_hinge" %in% run_sections) invisible(fit_temporal(tmsOcc_l, occ_formulas$prec_hinge, fit_tag("prec_hinge"), settings = fit_settings, inits = lambda_inits))

# Habitat x distance outside forest ----
if ("dist_hinge" %in% run_sections) invisible(fit_temporal(tmsOcc_l, occ_formulas$dist_hinge, fit_tag("dist_hinge"), settings = fit_settings, inits = lambda_inits))

# Habitat x time ----
if ("year_hinge" %in% run_sections) invisible(fit_temporal(tmsOcc_l, occ_formulas$year_hinge, fit_tag("year_hinge"), settings = fit_settings, inits = lambda_inits))

# Habitat x silvopasture cover within 300 m ----
if ("hab_ssp" %in% run_sections) invisible(fit_temporal(tmsOcc_l, occ_formulas$hab_ssp, fit_tag("hab_ssp"), settings = fit_settings, inits = lambda_inits))

# Silvopasture split by type ----
if ("split_hinge" %in% run_sections) invisible(fit_temporal(tmsOcc_l, occ_formulas$split_hinge, fit_tag("split_hinge"), settings = fit_settings, inits = lambda_inits))

# Habitat by landscape canopy cover ----
walk(intersect(c("ctx1500", "ctx1500_noht"), run_sections), \(model) fit_temporal(tmsOcc_l, occ_formulas[[model]], fit_tag(model), settings = fit_settings, inits = lambda_inits))

# Habitat-conversion models ----
walk(intersect(str_subset(names(occ_formulas), "^full1500"), run_sections), \(model) fit_temporal(tmsOcc_l, occ_formulas[[model]], fit_tag(model), settings = fit_settings, inits = lambda_inits))

# Additive model with silvopasture split by type ----
if ("add_split" %in% run_sections) invisible(fit_temporal(tmsOcc_l, occ_formulas$add_split, fit_tag("add_split"), settings = fit_settings, inits = lambda_inits))
