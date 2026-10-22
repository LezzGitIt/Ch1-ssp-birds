## PhD birds in silvopastoral landscapes##
# Analysis occupancy modeling -- Occupancy modeling using package spOccupancy

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(readxl)
library(sf)
library(unmarked)
library(spAbundance)
library(spOccupancy)
library(ubms) # function get_stancode() could be used to get code from model
library(flocker)
library(brms)
library(AICcmodavg)
library(furrr) # Use furrr::future_map()
library(coda)
library(broom.mixed)
library(tictoc)
library(ggpubr)
library(cowplot)
library(conflicted)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

## Data
#load("Rdata/Occ_abu_inputs_01.28.26.Rdata")
#load("Rdata/Occ_abu_models_01.19.25.Rdata")

# Create space
rm(list = ls()[!(ls() %in% c("msOcc_l", "ms_39_null", "ms_39_cc"))])

## Load data
Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"

## NOTE: Removing Hatico points until can digitize
Hatico_pc_nums <- str_pad(1:12, width = 2, pad = 0)

Pc_locs_dc_sf <- st_read(
  paste0(Wrangling_repo, "Derived_geospatial/shp/Pc_locs_dc.gpkg")
) %>% filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 

source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Rcookbook/Themes_funs.R")

# NOTES:
# 1) Requires unique spatial coordinates, so I added a small jitter (.001)
# 2) Does not allow NAs in occurrence (site) covs

# Define functions --------------------------------------------------------
# Extract parameter estimates and estimates of gof (ESS, rhat) from model list

extract_parms <- function(model, spatial = FALSE, ms = FALSE){
  
  beta <- tidyMCMC(model$beta.samples) %>% mutate(parameter = "Beta")
  alpha <- tidyMCMC(model$alpha.samples) %>% mutate(parameter = "Alpha")
  
  theta <- NULL
  if(spatial && !is.null(model$theta.samples)){
    theta <- tidyMCMC(model$theta.samples) %>% mutate(parameter = "Theta")
  }
  
  pars <- bind_rows(beta, alpha, theta)
  
  diag_tbl <- tibble(
    ess = c(model$ESS$beta, model$ESS$alpha, if(spatial) model$ESS$theta else NULL),
    rhat = c(model$rhat$beta, model$rhat$alpha, if(spatial) model$rhat$theta else NULL)
  )
  parms_df <- bind_cols(pars, diag_tbl)
  
  if(ms){
    parms_df <- parms_df %>% mutate(
      Species = str_split_i(term, "-", 2),
      term = str_split_i(term, "-", 1), 
      Species = ifelse(term == "phi", paste0("factor", Species), Species)
    ) %>% relocate(Species, .before = term)
  }
  return(parms_df)
}

# Rhat and effective sample size
plot_diagnostics <- function(parms_df){
  parms_df %>%
    ggplot(aes(ess, rhat, color = term, shape = parameter)) +
    geom_point() +
    geom_hline(yintercept = 1.1) +
    geom_vline(xintercept = 100) +
    guides(color = "none")
}

plot_species_diagnostics <- function(parms_df){
  parms_df %>%
    pivot_longer(cols = c(ess, rhat), names_to = "diagnostic") %>%
    ggplot(aes(species_ayerbe, value, color = term, shape = parameter)) +
    geom_point(position = "jitter") +
    facet_wrap(~diagnostic, scales = "free_y") +
    guides(color = "none") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

## Bayesian p-value
# Function taking posterior predictive check object and calculating bayesian p-value
calc_bayes_p <- function(ppc_obj){
  fit_stat <- ppc_obj$fit.stat
  group <- ifelse(ppc_obj$group == 1, "Site", "Replicate")
  
  tot.post <- ppc_obj$n.chains * ppc_obj$n.post
  bayes_p <- sum(ppc_obj$fit.y.rep > ppc_obj$fit.y) / tot.post
  bayes_p <- round(bayes_p, 3)
  tibble(fit_stat, group, bayes_p)
}

## GOF tests
# Run posterior predictive checks in spOccupancy
run_ppc <- function(mod){
  list(
    ppcOcc(mod, "freeman-tukey", 1),
    ppcOcc(mod, "freeman-tukey", 2),
    ppcOcc(mod, "chi-squared", 1),
    ppcOcc(mod, "chi-squared", 2)
  )
}

extract_gof <- function(ppc_list){
  imap(ppc_list, ~ map(.x, calc_bayes_p) %>%
         bind_rows() %>%
         mutate(species_ayerbe = .y)) %>%
    bind_rows()
}

# Visualize
plot_gof <- function(gof_tbl){
  gof_tbl %>%
    ggplot(aes(species_ayerbe, bayes_p, shape = group)) +
    geom_boxplot(position = position_dodge(.7)) +
    geom_point(aes(color = fit_stat), position = position_dodge(.7)) +
    geom_hline(yintercept = c(0.5, 0.1), linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "top")
}

plot_parm_estimates <- function(parms_df, ms = FALSE){
  parms_df <- parms_df %>%
    mutate(UCL = estimate + 1.96 * std.error,
           LCL = estimate - 1.96 * std.error)
  p <- ggplot(data = parms_df, aes(x = species_ayerbe, y = estimate, color = term))
  if(ms){ 
    p <- ggplot(data = parms_df, aes(x = Species, y = estimate, color = term))
  }
  p + geom_point(position = position_dodge(.7)) +
    geom_errorbar(position = position_dodge(.7), aes(ymin = LCL, ymax = UCL)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), 
          legend.position = "top") 
}

# Non-spatial model -------------------------------------------------------
safe_PGOcc <- safely(PGOcc)

# Random effect isn't the issue 
map_dbl(spOcc_ll, \(spOcc){
  spOcc$occ.covs %>% pull(Ecoregion) %>% 
    unique() %>% length()
}) %>% sort()
names(spOcc_ll)

Occ_mod_l_safe <- imap(spOcc_ll[-c(19)], \(spOcc, name){
  print(name)
  if(length(unique(spOcc$occ.covs$Ecoregion)) > 1){
    occ.formula <- as.formula("~ Ecoregion + Elev + Tot_prec + te + ssp + Habitat + Ano1 + (1 | Id_group_no_dc_num)") 
  } else{
    occ.formula <- as.formula("~ Elev + Tot_prec + te + ssp + Habitat + Ano1 + (1 | Id_group_no_dc_num)")
  }
  
  # Detection formula conditional on number of data collectors 
  if(length(unique(spOcc$det.covs$Nombre_institucion)) > 1){
    det.formula <- as.formula("~ Julian_day + Pc_start + log(Sampling_day)")
  } else{
    det.formula <- as.formula("~ Julian_day + Pc_start + log(Sampling_day)")
  }
  
  safe_PGOcc(
    occ.formula = occ.formula,
    det.formula = det.formula,
    data = spOcc[1:3], n.burn = 2000, n.samples = 4000,
    n.thin = 1, n.chains = 3, verbose = TRUE
    #k.fold = 2
  )
}, .progress = TRUE)

# Extract the successful mod runs
Occ_mod_l <- map(Occ_mod_l_safe, "result")
Occ_mod_l <- discard(Occ_mod_l, is.null) 

residuals(Occ_mod_l[[1]])

# Inspect errors
fit_errors <- map(Occ_mod_l_safe, "error")
discard(fit_errors, is.null) 

# >Extract parameters -----------------------------------------------------
parms_df <- extract_parms(Occ_mod_l, spatial = FALSE)

# >Diagnostics ----------------------------------------------------------------
# >>Rhat & ess ------------------------------------------------------------
# Visualize
plot_diagnostics(parms_df)
plot_species_diagnostics(parms_df)

# >>Goodness of fit ------------------------------------------------------
# Run custom functions to extract GOF 
gof_tbl <- Occ_mod_l %>% 
  run_ppc() %>% 
  extract_gof()

plot_gof(gof_tbl)

# >Inspect parm estimates -------------------------------------------------
parms_df %>% filter(parameter == "Beta", ess > 100, rhat < 1.1) %>% 
  plot_parm_estimates()

# Spatial model -----------------------------------------------------------
# >Priors -----------------------------------------------------------------
## Priors on spatial decay parameter (phi)
dist_mat <- st_distance(Pc_locs_dc_sf)
# By default the priors on phi are uniform(3 / max(dist_mat), 3 / min(dist_mat)), which is problematic as the upper boundary is undefined (division by 0). 
# Default priors
3 / max(dist_mat) # lower boundary
3 / min(dist_mat) # upper boundary

# Improved priors
lower <- 3 / 25e3 # largest distance between points within an ecoregion is 175 km
upper <- 3 / 500 # points are clustered, so replace min(dist_mat) with 'transect' length (500 m)
eff_spatial_range <- c(3 / upper, 3 / lower)
eff_spatial_range / 1000 # Translate above priors to the possible effective spatial ranges you are allowing (in km)

priors <- list(beta.normal = list(mean = 0, var = 2.72), 
               alpha.normal = list(mean = 0, var = 2.72), 
               sigma.sq.unif = c(0.5, 4), #c(.25, 2)
               phi.unif = c(lower, upper))

# error in base::stop("contrasts can be applied only to factors with 2 or more levels")
unique(spOcc_ll[[19]]$det.covs$Nombre_institucion)

# Create safe model fitting function to capture errors 
safe_spPGOcc <- safely(spPGOcc)

# spOcc_mod_l_safe_int, spOcc_mod_l_safe_re
spOcc_mod_l_safe_forest_1k <- imap(spOcc_ll, \(spOcc, spp_name){
  print(spp_name)
  
  # Occupancy formula conditional on number of ecoregions present
  #if(FALSE){
    if(length(unique(spOcc$occ.covs$Ecoregion)) > 1){
      occ.formula <- as.formula("~ Ecoregion + Elev + Tot_prec + te + ssp + Habitat + Canopy_cover + (1 | Id_group_no_dc_num)") 
    } else{
      occ.formula <- as.formula("~ Elev + Tot_prec + te + ssp + Habitat + Canopy_cover + (1 | Id_group_no_dc_num)")
    }
  
  #if(FALSE){
    # Detection formula conditional on number of data collectors 
    if(length(unique(spOcc$det.covs$Nombre_institucion)) > 1){
      det.formula <- as.formula("~ Nombre_institucion + Julian_day + Pc_start")
    } else{
      det.formula <- as.formula("~ Julian_day + Pc_start")
    }
  #}
  
  safe_spPGOcc(
    occ.formula = occ.formula, 
    det.formula = det.formula,
    priors = priors,
    cov.model = "exponential", NNGP = TRUE, n.neighbors = 15,
    data = spOcc[1:4], n.burn = 3000, n.batch = 400, batch.length = 25,
    n.thin = 20, n.chains = 3, n.report = 100, verbose = TRUE
  )
}, .progress = TRUE)

# Extract results from successful runs
spOcc_mod_l_forest_1k <- map(spOcc_mod_l_safe_forest_1k, "result") # unif still has 3 / .1 for upper
spOcc_mod_l_forest_1k <- discard(spOcc_mod_l_forest_1k, is.null)
# Inspect errors
fit_errors <- map(spOcc_mod_l_safe_forest_1k, "error")
discard(fit_errors, is.null) 

# >Extract parameters -----------------------------------------------------
mod_l <- spOcc_mod_l_forest_1k #int, re, eco_re, eco

parms_df <- map(mod_l, extract_parms, spatial = TRUE)

# >Diagnostics ----------------------------------------------------------------
# >>Rhat & ess ------------------------------------------------------------
# Visualize
parms_df %>% filter(ess < 1000 & rhat < 10) %>%
  plot_diagnostics()
parms_df %>% filter(ess < 1000) %>%
  plot_species_diagnostics()

# >>Goodness of fit ------------------------------------------------------
# Run custom functions to extract GOF 
# These ppc models don't fit unless you use n.thin
gof_tbl_forest_1k <- spOcc_mod_l_forest_1k %>%  
  run_ppc() %>% 
  extract_gof()

plot_gof(gof_tbl_forest_1k)

# Traceplots
imap(spOcc_mod_l_forest_1k, \(mod, name){
  plot(mod, param = "theta")
})

# >Inspect parm estimates -------------------------------------------------
parms_df %>% 
  filter(parameter == "Theta" & term == "phi") %>% # & term == "phi" sigma.sq
  plot_parm_estimates()

# Occupancy covariates
imap(spOcc_mod_l_forest_1k, \(mod, name){
  MCMCvis::MCMCplot(
    mod$beta.samples, rank = TRUE, main = name, ref_ovl = TRUE, exact = FALSE, excl = "Ecoregion", params = "Habitat" # "all"
    )
})

# Detection covariates
imap(spOcc_mod_l_forest_1k, \(mod, name){
  MCMCvis::MCMCplot(mod$alpha.samples, rank = TRUE, main = name, ref_ovl = TRUE)
})

# Can play with bayesplot::mcmc_dens() function as well, or as.data.frame() to get posterior draws 

# >Effective spatial range ------------------------------------------------
# Examine the estimated effective spatial ranges
effective_sp_range <- parms_df$Cyanocorax_violaceus %>% 
  filter(parameter == "Theta" & term == "sigma.sq") %>% #"sigma.sq" "phi"
  rename(phi_estimate = estimate) %>%
  select(-c(parameter, term)) %>%
  mutate(eff_range_km = (3 / phi_estimate) / 1000)
effective_sp_range %>% pull(eff_range_km) %>% 
  sort() %>%
  setNames(NULL)

# Look for relationship between fit and effective range
effective_sp_range %>% #filter(rhat < 100) %>%
  ggplot(aes(x = ess, y = eff_range_km)) +
  geom_point() + 
  geom_smooth()

# Model selection ---------------------------------------------------------
## Compare spatial vs non-spatial
# Put model lists into a named list
mods_comp_ll <- list(spOcc_mod_l_forest_1k)
mods_comp_ll <- setNames(mods_comp_ll, c("forest_1k"))

# Nested map call to create model comparison tibble 
mods_comp <- map(mods_comp_ll, \(mod_l){
  map(mod_l, ~ waicOcc(.x)) %>% bind_rows(.id = "Species_ayerbe_")
}) %>% bind_rows(.id = "Model")

# Visualize
mods_comp %>% ggplot(aes(x = Model, y = WAIC)) +
  geom_boxplot() +
  geom_point(aes(color = Species_ayerbe_)) + 
  guides(color = "none")

# Examine mean WAIC scores across species
mods_comp %>% summarize(mean = mean(WAIC), .by = Model) %>% 
  arrange(mean)

# EXTRAS ------------------------------------------------------------------

# Plot ppc
ppc.df <- data.frame(fit = ppc_sites$fit.y,
                     fit.rep = ppc_sites$fit.y.rep,
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, 
     pch = 21, ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')

diff.fit <- ppc_sites$fit.y.rep.group.quants[3, ] - ppc_sites$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')

names(spOcc_mod[[1]])
summary(spOcc_mod[[3]])
map(spOcc_mod, summary)

# Multi-season models
?stPGOcc # Single species
?stMsPGOcc # Multi species
?svcTMsPGOcc # Svc multi species

## We want a model that is... In the following list, numbers 1 & 2 are our inferential goals, 3, 4, & 5 are not primary inference goals but are likely important to accurately estimate the effects we are interested in. 
# 1) svc: Allows for parameter estimates to vary through space, good as we want to nail the parameter estimates for silvopasture, but have 5 distinct (& heterogeneous) ecoregions. From vignette: 'The model reduces to a traditional single-species occupancy model when all covariate eﬀects are assumed constant across space and a spatial occupancy model (Johnson et al. 2013; Doser et al. 2022) when only the intercept is assumed to vary across space'. So this is essentially random slope & random intercept
# 2) multi-season: Important to understand how species occupancy varies through time in the different habitat types (forest, silvopasture, pasture)
# 3) spatial: Account for spatial autocorrelation. From vignette, this is essentially a random intercept
# 4) latent factor: Account for residual species correlations via a factor modeling approach, which is a dimension reduction technique that can account for correlations among a large number of species without a massive computational cost
# 5) integrated: We have different data collectors with different abilities to detect birds by sight & (especially) by ear

# NOTE there are several other relevant arguments if we add coords back into the mix (see below)

# spAbundance ------------------------------------------------------------
# Doser suggests that including spatial autocorrelation & negative-binomial distribution will help with zero-inflation

# Useful functions
?spAbundance::sfMsNMix
?NMix
postHocLM # Explore how species traits relate to species-specific parameter estimates from a multi-species occupancy model

# Save --------------------------------------------------------------------
# Save model output
#save(fit_ml, fit_ml_safe, fit_stan, Error_l, Error_meta, Ls_count, file = paste0("Rdata/Occ_abu_models_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))