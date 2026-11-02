## PhD birds in silvopastoral landscapes ##
# Multi-season (temporal) spatial multi-species occupancy model (spOccupancy::svcTMsPGOcc with svc.cols = 1, i.e. a spatially varying intercept only), adapted from Doser's space-time vignette (Extra_scripts/spaceTimeModelsHTML.Rmd): each physical location is one site with one spatial effect, and the year-group x season survey periods are a separate time dimension, so stacked copies of a location no longer need jittered coordinates
# Model family: Bernoulli occupancy and Bernoulli detection with logit links (Polya-Gamma sampler); spatial NNGP exponential process with 6 latent factors; no AR(1) temporal effect (periods are unequally spaced 2013-2026), the trend is carried by the Year1 fixed effect as in the stacked model
# Non-spatial comparison (spOccupancy::tMsPGOcc): same likelihood and links, no spatial process and no latent factors (tMsPGOcc has neither); periods at a location are linked only through an unstructured location random effect, as in the vignette's tPGOcc example
# svcTMsPGOcc with svc.cols = 1 is the same model as stMsPGOcc (seeded runs give bit-identical draws, 2026-11-02) but, unlike stMsPGOcc, accepts range.ind (species x location, so clipping is spatial, not seasonal); tMsPGOcc has no range.ind, so the non-spatial fits are unclipped
# Run headless with: nohup caffeinate -is Rscript Scripts/spOccupancy_temporal.R <section(s)> > <log> 2>&1 &  (sections: additive, prec, prec_nofarm, prec_nofarm_long, prec_ns, prec_ns_loc, year_nofarm_long, dist_nofarm_long; or smoke for a 10-batch, 1-chain dimension check)

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(spOccupancy)
# Pin functions that other attached packages mask (a plain assignment rather than conflicted, which turns every ambiguous call into an error)
select <- dplyr::select
filter <- dplyr::filter

# Configuration -----------------------------------------------------------
## Stacked list the temporal list is built from (one row per location x year-group x season)
msOcc_l_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_278_hab_fix.rds"
tmsOcc_l_path <- "Derived/Rdata/Model_input/Temporal/tmsOcc_l_meta_278_hab_fix.rds"

## Site = physical location; time = year-group x season in chronological order (replicates stay the within-season visits, so closure is the same as in the stacked model)
Site_col <- "Id_survey_no_dc"
Time_cols <- c("Year_grp", "Season")
Time_levels <- c("13-14.Early", "13-14.Late", "16-17.Early", "16-17.Late", "19.Early", "22.Early", "24.Early", "25-26.Early")
## The multi-season models need occupancy covariates in every site x period cell, including periods a location was not surveyed (no detection data there): covariates that describe the period are filled with the period mean, all other time-varying covariates with the location's value from its nearest surveyed period
Period_covs <- "Year1"

## Sections to run; command-line arguments override this
run_sections <- c("additive", "prec")
if (length(commandArgs(trailingOnly = TRUE)) > 0) run_sections <- commandArgs(trailingOnly = TRUE)

## Formulas
occ_formulas <- list(
  additive = ~ Habitat + ssp + Dist_forest + Year1 + Tot_prec + (1 | Farm_join_num),
  prec = ~ Habitat * Tot_prec + Dist_forest + ssp + Year1 + (1 | Farm_join_num),
  # Spatial process and Farm_join compete (intercept and farm variance mix poorly together, model log rows 71-72 and 77-78), so each is fitted without the other
  prec_nofarm = ~ Habitat * Tot_prec + Dist_forest + ssp + Year1,
  prec_ns = ~ Habitat * Tot_prec + Dist_forest + ssp + Year1 + (1 | Farm_join_num) + (1 | Site_num),
  # Location and Farm_join variances compete in prec_ns (locations are nested in the 13 units), so the location effect is also fitted alone
  prec_ns_loc = ~ Habitat * Tot_prec + Dist_forest + ssp + Year1 + (1 | Site_num),
  # One habitat interaction per model (two in one model judged too costly), each on the primary no-Farm_join structure
  year_nofarm = ~ Habitat * Year1 + Tot_prec + Dist_forest + ssp,
  dist_nofarm = ~ Habitat * Dist_forest + Tot_prec + ssp + Year1
)
det.formula <- ~ Forest_bin + Pc_start + I(Pc_start^2) + Prob_breed + Sampling_day_log + (1 | Collector_team_num)

## MCMC settings (full length, as every logged stacked fit; the screening length did not predict mixing problems)
# Thinned by 25, not 10 as in the stacked fits: the species x site x period posterior arrays (z, psi, likelihood) are 3.6x the stacked size, and at 3,000 kept samples assembling them exceeded R's 18 GB memory limit after sampling (2026-10-30); 1,200 kept samples (400 per chain; n.thin must divide the 10,000 post-burn samples) keep them near stacked-fit size
batch.length <- 25
mcmc_full <- list(n.batch = 1200, n.burn = 20000, n.thin = 25, n.chains = 3)
# Twice the post-burn-in length for parameters still just above Rhat 1.1; thinned by 75 so the kept samples (400 per chain) and memory stay as in mcmc_full
mcmc_long <- list(n.batch = 2400, n.burn = 30000, n.thin = 75, n.chains = 3)
n.factors <- 6; n.neighbors <- 5
priors <- list(phi.unif = list(3 / 25e3, 3 / 23e3)) # effective spatial range 23-25 km, as in the stacked fits
ar1 <- FALSE

# Functions ---------------------------------------------------------------
# Fit one multi-season model and save it under a tag; each call is independent so a crashed or cached fit never affects another
fit_temporal <- function(data, occ.formula, tag, spatial = TRUE, mcmc = mcmc_full, smoke = FALSE) {
  # Latent occurrence starting values from the vignette: 1 where the species was detected in that site x period, else 0
  z.inits <- apply(data$y, c(1, 2, 3), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
  len <- if (smoke) list(n.batch = 10, n.burn = 0, n.thin = 1, n.chains = 1) else mcmc
  save_path <- paste0("Derived/Rdata/Model_output/Temporal/Occ_mod_", dim(data$y)[1], "spp_1ecor_temporal_", tag, "_", Sys.Date(), ".rds")
  if (!smoke && file.exists(save_path)) stop("'", save_path, "' already exists")
  # Create the output folder up front so a long fit never fails at the save step on a fresh machine
  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  cat("Started", tag, format(Sys.time()), "\n")
  start <- Sys.time()
  mod <- if (spatial) {
    svcTMsPGOcc(occ.formula = occ.formula, det.formula = det.formula, data = data, inits = list(z = z.inits), priors = priors, tuning = list(phi = 2), svc.cols = 1, cov.model = "exponential", NNGP = TRUE, n.neighbors = n.neighbors, n.factors = n.factors, n.batch = len$n.batch, batch.length = batch.length, accept.rate = 0.43, ar1 = ar1, n.omp.threads = 1, verbose = TRUE, n.report = 100, n.burn = len$n.burn, n.thin = len$n.thin, n.chains = len$n.chains)
  } else {
    tMsPGOcc(occ.formula = occ.formula, det.formula = det.formula, data = data[c("y", "occ.covs", "det.covs")], inits = list(z = z.inits), n.batch = len$n.batch, batch.length = batch.length, accept.rate = 0.43, ar1 = ar1, n.omp.threads = 1, verbose = TRUE, n.report = 100, n.burn = len$n.burn, n.thin = len$n.thin, n.chains = len$n.chains)
  }
  cat("Finished", tag, format(Sys.time()), "| wall", format(Sys.time() - start), "\n")
  if (!smoke) { saveRDS(mod, save_path); cat("SAVED", save_path, "\n") }
  mod
}

# Build the temporal list -------------------------------------------------
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
      if (nm %in% Period_covs) return(matrix(colMeans(m, na.rm = TRUE), n_site, n_time, byrow = TRUE, dimnames = dimnames(m)))
      t(apply(m, 1, function(r) { obs <- which(!is.na(r)); r[sapply(seq_along(r), function(t) obs[which.min(abs(obs - t))])] }))
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
## Location index for the non-spatial model's location random effect (not stored in the list; unused by the spatial model)
tmsOcc_l$occ.covs$Site_num <- seq_len(dim(tmsOcc_l$y)[2])
str(tmsOcc_l, max.level = 2)

# Smoke test --------------------------------------------------------------
# Catches formula / dimension errors in about a minute; not for judging convergence
if ("smoke" %in% run_sections) smoke_mod <- fit_temporal(tmsOcc_l, occ_formulas$prec, "smoke", smoke = TRUE)

# Additive ----------------------------------------------------------------
# Same formula and data as the converged stacked fit D (Occ_mod_98spp_1ecor_add_ssp_habfix_bfcoords_2026-10-30.rds), so WAIC and estimates compare directly
if ("additive" %in% run_sections) Occ_mod_add <- fit_temporal(tmsOcc_l, occ_formulas$additive, "add")

# Habitat x Tot_prec ------------------------------------------------------
# Primary formula; the stacked version mixed poorly at full length (model log row 70)
if ("prec" %in% run_sections) Occ_mod_prec <- fit_temporal(tmsOcc_l, occ_formulas$prec, "prec")

# Habitat x Tot_prec, spatial, no Farm_join ------------------------------
# Tests whether dropping the farm effect lets the spatial model mix (as it did for the stacked grid.index fit, model log row 78)
if ("prec_nofarm" %in% run_sections) Occ_mod_prec_nofarm <- fit_temporal(tmsOcc_l, occ_formulas$prec_nofarm, "prec_nofarm")

# Habitat x Tot_prec, non-spatial, Farm_join + location ---------------------
# The other half of the test: farm effect kept, spatial process replaced by an unstructured location effect
if ("prec_ns" %in% run_sections) Occ_mod_prec_ns <- fit_temporal(tmsOcc_l, occ_formulas$prec_ns, "prec_ns", spatial = FALSE)

# Habitat x Tot_prec, spatial, no Farm_join, long run ---------------------
# Same model as prec_nofarm at twice the post-burn-in length, to see whether Tot_prec (Rhat 1.14, model log row 83) converges with more samples
if ("prec_nofarm_long" %in% run_sections) Occ_mod_prec_nofarm_long <- fit_temporal(tmsOcc_l, occ_formulas$prec_nofarm, "prec_nofarm_long", mcmc = mcmc_long)

# Habitat x Tot_prec, non-spatial, location only --------------------------
if ("prec_ns_loc" %in% run_sections) Occ_mod_prec_ns_loc <- fit_temporal(tmsOcc_l, occ_formulas$prec_ns_loc, "prec_ns_loc", spatial = FALSE)

# Habitat x Year1, spatial, no Farm_join, long run ------------------------
if ("year_nofarm_long" %in% run_sections) Occ_mod_year_nofarm_long <- fit_temporal(tmsOcc_l, occ_formulas$year_nofarm, "year_nofarm_long", mcmc = mcmc_long)

# Habitat x Dist_forest, spatial, no Farm_join, long run ------------------
if ("dist_nofarm_long" %in% run_sections) Occ_mod_dist_nofarm_long <- fit_temporal(tmsOcc_l, occ_formulas$dist_nofarm, "dist_nofarm_long", mcmc = mcmc_long)
