## PhD birds in silvopastoral landscapes ##
# Multi-season (temporal) spatial multi-species occupancy model (spOccupancy::svcTMsPGOcc with svc.cols = 1, i.e. a spatially varying intercept only), adapted from Doser's space-time vignette (Extra_scripts/spaceTimeModelsHTML.Rmd): each physical location is one site with one spatial effect, and the year-group x season survey periods are a separate time dimension, so stacked copies of a location no longer need jittered coordinates
# Model family: Bernoulli occupancy and Bernoulli detection with logit links (Polya-Gamma sampler); spatial NNGP exponential process with 2 latent factors, their loadings started at an earlier fit's posterior means; no AR(1) temporal effect (periods are unequally spaced 2013-2026), time enters as a fixed effect
# Distance to the forest edge is hinge-coded: Dist_open (distance outside forest, 0 in forest) and Dist_inside (depth inside forest, 0 outside), so each side of the edge has its own slope and the habitat contrasts sit at the edge
# svcTMsPGOcc with svc.cols = 1 is the same model as stMsPGOcc (seeded runs give bit-identical draws, 2026-11-02) but, unlike stMsPGOcc, accepts range.ind (species x location, so clipping is spatial, not seasonal)
# Run headless with: nohup caffeinate -is Rscript Scripts/spOccupancy_temporal.R <section(s)> > <log> 2>&1 &  (sections: prec_hinge, dist_hinge, year_hinge, hab_ssp, split_hinge; or smoke for a 10-batch, 1-chain dimension check of every model)

# Load libraries ----
library(tidyverse)
library(spOccupancy)
# Pin functions that other attached packages mask (a plain assignment rather than conflicted, which turns every ambiguous call into an error)
select <- dplyr::select
filter <- dplyr::filter
source("Scripts/spOcc_fns.R")

# Configuration ----
## Stacked list the temporal list is built from (one row per location x year-group x season), with the scaling of its covariates
msOcc_l_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm.rds"
tmsOcc_l_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_98_hab_fix.rds"
Scaling_path <- str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv")
## Farm_diversity's landscape metrics (forest = % forest within 300 m, from the land cover matched to the survey year; survey ID to location code), the wrangling repo's site table (silvopasture type per location) and visit table (survey dates)
Lsm_path <- "../Farm_diversity/Derived/Excels/Event_covs_lsm.csv"
Site_covs_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Site_covs.csv"
Event_covs_path <- "../Ssp-bird-data-wrangling/Derived/Excels/Event_covs.csv"

## Site = physical location; time = year-group x season in chronological order (replicates stay the within-season visits, so closure is the same as in the stacked model)
Site_col <- "Id_survey_no_dc"
Time_cols <- c("Year_grp", "Season")
Time_levels <- c("13-14.Early", "13-14.Late", "16-17.Early", "16-17.Late", "19.Early", "22.Early", "24.Early", "25-26.Early")
## The multi-season models need occupancy covariates in every site x period cell, including periods a location was not surveyed (no detection data there): covariates that describe the period are filled with the period mean, all other time-varying covariates with the location's value from its nearest surveyed period
Period_covs <- "Survey_date"

## Sections to run, from the command line; none when the script is sourced (builds the data only)
run_sections <- commandArgs(trailingOnly = TRUE)

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
fit_tag <- \(model) paste0(model, if (Time_cov == "Survey_date") "_date", "_f2_inits")
det.formula <- ~ Forest_bin + Pc_start + I(Pc_start^2) + Prob_breed + Sampling_day_log + (1 | Collector_team_num)

## MCMC settings: 60,000 iterations, 30,000 burn-in, thinned by 75 (400 kept per chain); the species x site x period posterior arrays (z, psi, likelihood) are 3.6x the stacked size, and at 3,000 kept samples assembling them exceeded R's 18 GB memory limit (2026-10-30)
batch.length <- 25
mcmc_long <- list(n.batch = 2400, n.burn = 30000, n.thin = 75, n.chains = 3)
n.neighbors <- 5
priors <- list(phi.unif = list(3 / 25e3, 3 / 23e3)) # effective spatial range 23-25 km
ar1 <- FALSE
## 2 spatial factors: with 3, factor 2 emptied even when started informative (Multi_spp_100 row 110); started at the posterior-mean loadings of the 2-factor signed-distance fit (row 90), where factor 2 was informative, because 2-factor runs started at 0 often collapsed it (Jeff Doser's remedy for loadings that won't converge); spOccupancy uses the same loadings start in every chain
Lambda_source_path <- "Derived/Rdata/Model_output/Temporal/Occ_mod_98spp_prec_nofarm_long_f2_2026-11-05.rds"
Lambda_inits_path <- "Derived/Rdata/Model_output/Temporal/lambda_inits_prec_nofarm_long_f2.rds"
## Settings passed to fit_temporal() (spOcc_fns.R); every fit is seeded so a rerun reproduces its draws exactly
fit_settings <- list(det.formula = det.formula, priors = priors, n.neighbors = n.neighbors, n.factors = 2, batch.length = batch.length, ar1 = ar1, mcmc = mcmc_long, out_dir = "Derived/Rdata/Model_output/Temporal", seed = 20261107)

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
tmsOcc_l$occ.covs <- oc
str(tmsOcc_l, max.level = 2)

# Loadings starting values ----
## Posterior means of the earlier fit's loadings (species x factor), computed once; matched to species by position, so the species order must equal the data's
if (!file.exists(Lambda_inits_path)) {
  Lambda_source <- readRDS(Lambda_source_path)
  saveRDS(matrix(colMeans(Lambda_source$lambda.samples), nrow = length(Lambda_source$sp.names), dimnames = list(Lambda_source$sp.names, NULL)), Lambda_inits_path)
  rm(Lambda_source); gc()
}
lambda_inits <- list(lambda = list(readRDS(Lambda_inits_path)))
stopifnot(identical(rownames(lambda_inits$lambda[[1]]), dimnames(tmsOcc_l$y)[[1]]))

# Smoke test ----
# Catches formula / dimension errors in about a minute per model; not for judging convergence
if ("smoke" %in% run_sections) iwalk(occ_formulas, \(f, model) fit_temporal(tmsOcc_l, f, fit_tag(model), settings = fit_settings, inits = lambda_inits, smoke = TRUE))

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
