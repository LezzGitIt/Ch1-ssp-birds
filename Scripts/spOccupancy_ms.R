## PhD birds in silvopastoral landscapes##
# Analysis multi-species occupancy modeling -- Occupancy modeling using package spOccupancy

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(readxl)
library(sf)
library(spOccupancy)
library(coda)
library(broom.mixed)
library(ggpubr)
library(cowplot)
library(conflicted)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

## Data
#load("Rdata/Occ_abu_inputs_01.28.26.Rdata")
#load("Rdata/Occ_abu_models_01.19.25.Rdata")

# Create space by removing old models
rm(list = ls()[!(ls() %in% c())])
msOcc_l <- readRDS("Rdata/msOcc_l_season.rds")

## Load data
Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"

## NOTE: Removing Hatico points until can digitize
Hatico_pc_nums <- str_pad(1:12, width = 2, pad = 0)

Pc_locs_dc_sf <- st_read(
  paste0(Wrangling_repo, "Derived_geospatial/shp/Pc_locs_dc.gpkg")
) %>% filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 

source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")
source("Scripts/spOcc_fns.R")

# NOTES:
# 1) Requires unique spatial coordinates, so I added a small jitter (.001)
# 2) Does not allow NAs in occurrence (site) covs

# Row identifier ----------------------------------------------------------
#NOTE:: The unique row identifier  [e.g. ID x Year] is critical , this defines what a row is in your dataframes and how many rows each dataframe will have 
Row_identifier <- c("Id_muestreo", "Ano_grp", "Season") #, "Season"

# Non-spatial -----------------------------------------------------------
occ.formula <- as.formula("~ Ecoregion + Habitat + Ano1 + Canopy_cover + Tot_prec + te + ssp + (1 | Id_group_no_dc)") 
det.formula <- as.formula("~ Julian_day + Pc_start + Nombre_institucion + Forest_bin") 

# >Run mod ----------------------------------------------------------------
stop()
start <- Sys.time()
ms_160_ns <- lfMsPGOcc(
  occ.formula = occ.formula, 
  det.formula = det.formula, 
  data = msOcc_l, 
  n.samples = 5000,
  n.factors = 7,
  #inits = inits, 
  #priors = priors, 
  #tuning = tuning, 
  verbose = TRUE, 
  n.report = 100, 
  n.burn = 3000, 
  n.thin = 2, 
  n.chains = 3
)
end <- Sys.time()
end - start
Occ_mod <- ms_117_ns

# Expected log pointwise predictive density (elpd), the effective number of parameters (pD), waic is the sum of the waic values for each species 
waic_spp <- waicOcc(Occ_mod) 
waic_spp

# >Extract parameters -----------------------------------------------------
parms_df <- Occ_mod %>% 
  extract_parms(spatial = TRUE, ms = TRUE) 

# >Diagnostics ----------------------------------------------------------------
# >>Rhat & ess ------------------------------------------------------------
# Visualize
parms_df %>% plot_diagnostics()

# >>Goodness of fit ------------------------------------------------------
## Need to reduce the number of posterior samples for GOF check

# Logical true false vector for subsetting
TF <- str_detect(names(ms_100_8fac), "samples")
Samples_l <- ms_100_8fac[TF][-15]

# Subset 500 samples irrespective of dimensionality of the samples object
samples_500 <- map(Samples_l, \(samples) {
  random_sample <- sample(dim(samples)[1], 500, replace = FALSE)
  if(length(dim(samples)) == 2){
    samples[random_sample, , drop = FALSE]
  } 
  else if(length(dim(samples)) == 3){
    samples[random_sample, , , drop = FALSE]
  }
  else if(length(dim(samples)) == 4){
    samples[random_sample, , , , drop = FALSE]
  }
})

# Create new object 
samples_500_l <- ms_100_8fac
# Overwrite with the smaller subset of samples
samples_500_l[TF][-15] <- samples_500
samples_500_l$n.post <- 500
#samples_500_l$n.samples <- 500
class(samples_500_l) # class = 'svcMsPGOcc'

#rm(list = ls()[!(ls() %in% c("parms_df", "samples_500_l", "run_ppc"))])

# All four combinations produce the same error
ms_100_ppc <- ppcOcc(samples_500_l, fit.stat = "chi-squared", group = 1)

# Run custom functions to extract GOF 
# These ppc don't fit unless you use n.thin
ms_100_ppc <- ms_39_ls_hab %>%  
  run_ppc() #%>% 
gof_tbl_39 <- ms_100_ppc %>% extract_gof()

plot_gof(gof_tbl_39)

# Traceplots
plot(Occ_mod, param = "theta")

# >Inspect parm estimates -------------------------------------------------
parms_df %>% 
  filter(parameter == "Theta" & term == "phi") %>% # & term == "phi" sigma.sq
  plot_parm_estimates(ms = TRUE)

parms_df %>% filter(parameter == "Beta") %>% 
  filter(str_detect(term, "Habitat")) %>%
  plot_parm_estimates(ms = TRUE)

library(posterior)
library(ggdist)
draws_beta <- as_draws_df(Occ_mod$beta.comm.samples) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "value")
draws_alpha <- as_draws_df(Occ_mod$alpha.comm.samples) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "value")


# Only 1. chain? 
draws_beta %>% filter(param == ".chain") %>% 
  tabyl(value)

plot_half_eye <- function(df){
  df %>% 
    filter(!str_detect(param, "\\.")) %>% # rm .iteration, .draw, and .chain
    ggplot(aes(x = value, y = param, fill = param)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(alpha = .8) +
    labs(x = NULL, y = NULL) +
    guides(fill = "none")
}

# Understand what is contained in (Intercept)
draws_beta %>% 
  #filter(str_detect(param, "Habitat")) %>%
  filter(!str_detect(param, "Ecoregion|Intercept")) %>% 
  plot_half_eye() + 
  labs(title = "117 species")
draws_alpha %>% 
  filter(!str_detect(param, "institucion|Intercept")) %>%
  plot_half_eye() + 
  labs(title = "117 species")

ggsave("Figures/alpha_post_8fac.png", bg = "white") #beta

# Can play with bayesplot::mcmc_dens() function as well, or as.data.frame() to get posterior draws 

# Spatial -----------------------------------------------------------
occ.formula <- as.formula("~ Ecoregion + Habitat + Ano1 + Canopy_cover + Tot_prec + te + ssp + (1 | Id_group_no_dc)") 
det.formula <- as.formula("~ Forest_bin + Pc_start + Nombre_institucion + Prob_breed") 

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
               sigma.sq.unif = c(.5, 4), #c(.25, 2)
               phi.unif = list(lower, upper))


## Initial values
n.factors <- 7
# Number of species
N <- nrow(msOcc_l$y)
# Initiate all lambda initial values to 0. 
lambda.inits <- matrix(0, N, n.factors)
# Set diagonal elements to 1
diag(lambda.inits) <- 1
# Set lower triangular elements to random values from a standard normal distribution
N_to_fill <- sum(lower.tri(lambda.inits))
lambda.inits[lower.tri(lambda.inits)] <- rnorm(n = N_to_fill, mean = 0, sd = 1)
# Default initial values for lambda.
lambda.inits

## Run multi-species model accounting for spatial autocorrelation, detection, and biogeographic clipping
# The function svcMsPGOcc takes range.ind as an argument. range.ind is a matrix with rows corresponding to species and columns corresponding to sites, with each element taking value 1 if that site is within the range of the corresponding species and 0 if it is outside of the range.

# >Run mod ----------------------------------------------------------------
## Model paramaters
n.batch <- 200
batch.length <- 25
n.chains <- 3
Samples_per_chain <- n.batch * batch.length # Total number of samples 

# We throw away much of the MCMC work 
n.burn <- 3000 
n.thin <- 4

# Total samples per chain - burn in per chain
Keep_per_chain <- (Samples_per_chain - n.burn) / n.thin 
# This is how many samples we end up with in the end
Tot_samples <- Keep_per_chain * n.chains
Tot_samples

# dim(msOcc_l$y)
start <- Sys.time()
ms_100_season <- svcMsPGOcc(
  occ.formula = occ.formula, 
  det.formula = det.formula, 
  data = msOcc_l, 
  #inits = inits, 
  n.batch = n.batch, 
  batch.length = batch.length, 
  accept.rate = 0.43, 
  priors = priors, 
  n.factors = n.factors, # n.factors,
  svc.cols = 1,
  # model changes the tuning values after each batch of the MCMC to yield acceptance rates that are close to our target acceptance rate
  tuning = list(phi = 2), 
  n.omp.threads = 1, 
  verbose = TRUE, 
  n.report = 100, 
  n.burn = n.burn, 
  n.thin = n.thin, 
  n.chains = n.chains,
  # spatial autocorrelation
  NNGP = TRUE, 
  n.neighbors = 5, 
  cov.model = "exponential"
)
end <- Sys.time()
end - start
Occ_mod <- ms_100_season

# Expected log pointwise predictive density (elpd), the effective number of parameters (pD), waic is the sum of the waic values for each species 
waic_spp <- waicOcc(Occ_mod) 
waic_spp
stop()

# >Extract parameters -----------------------------------------------------
parms_df <- Occ_mod %>% 
  extract_parms(spatial = TRUE, ms = TRUE) %>% 
  left_join(Spp_ordered, by = join_by(Species_num == Order))

# >Diagnostics ----------------------------------------------------------------
summary(Occ_mod)

# >>Rhat & ess ------------------------------------------------------------
# Visualize
parms_df %>% plot_diagnostics()

# Which parameters are failing to converge? 
Prob_parms <- parms_df %>% 
  filter(rhat > 1.1 | ess < 100) %>% 
  arrange(desc(rhat))
Prob_parms %>% filter(str_detect(term, "Ecoregion")) 
Prob_parms %>% tabyl(parameter, term)

# Proportion of parameters that are problematic
round(nrow(Prob_parms) / nrow(parms_df), 2) # 12%

# >>Goodness of fit ------------------------------------------------------
## Need to reduce the number of posterior samples for GOF check

# Logical true false vector for subsetting
TF <- str_detect(names(ms_100_8fac), "samples")
Samples_l <- ms_100_8fac[TF][-15]

# Subset 500 samples irrespective of dimensionality of the samples object
samples_500 <- map(Samples_l, \(samples) {
  random_sample <- sample(dim(samples)[1], 500, replace = FALSE)
  if(length(dim(samples)) == 2){
    samples[random_sample, , drop = FALSE]
  } 
  else if(length(dim(samples)) == 3){
    samples[random_sample, , , drop = FALSE]
  }
  else if(length(dim(samples)) == 4){
    samples[random_sample, , , , drop = FALSE]
  }
})

# Create new object 
samples_500_l <- ms_100_8fac
# Overwrite with the smaller subset of samples
samples_500_l[TF][-15] <- samples_500
samples_500_l$n.post <- 500
#samples_500_l$n.samples <- 500
class(samples_500_l) # class = 'svcMsPGOcc'

#rm(list = ls()[!(ls() %in% c("parms_df", "samples_500_l", "run_ppc"))])

# All four combinations produce the same error
ms_100_ppc <- ppcOcc(samples_500_l, fit.stat = "chi-squared", group = 1)

# Run custom functions to extract GOF 
# These ppc don't fit unless you use n.thin
ms_100_ppc <- ms_39_ls_hab %>%  
  run_ppc() #%>% 
gof_tbl_39 <- ms_100_ppc %>% extract_gof()

plot_gof(gof_tbl_39)

# Traceplots - spatial parameters
plot(Occ_mod, param = "theta", density = FALSE)
# Traceplots - factor loadings 
#plot(Occ_mod$lambda.samples[,140:150], param = "lambda", density = FALSE)
plot(Occ_mod, param = "lambda", density = FALSE)

# >>Inspect factor loadings ------------------------------------------------
# Lambda samples are an mcmc object with dimensions of number of samples x (spp * number of factors). 
# The flat green lines are when the intercept is fixed at zero? 
Occ_mod$lambda.samples %>% dim()

sp_fac_loadings <- colSums(Occ_mod$lambda.samples)
# Extract species and latent factor number
species <- str_split_i(names(sp_fac_loadings), "-", 1)
factor <- str_split_i(names(sp_fac_loadings), "-", 2)
sp_fac_tbl <- tibble(species, factor, value = sp_fac_loadings)

# Doser suggests using the species factor loadings to set order 
sp_fac_tbl %>%
  ggplot() + geom_density(aes(x = value, color = factor))
sp_fac_tbl %>% mutate(abs_val = abs(value)) %>% 
  slice_max(abs_val, by = factor) 

# >Inspect parm estimates -------------------------------------------------
# Species-specific effects 
parms_df %>% 
  filter(parameter == "Theta" & term == "phi") %>% # & term == "phi" sigma.sq
  plot_parm_estimates(ms = TRUE)

parms_df %>% filter(parameter == "Beta") %>% 
  filter(str_detect(term, "Habitat")) %>%
  plot_parm_estimates(ms = TRUE)

parms_df %>% filter(parameter == "Alpha") %>% 
  filter(str_detect(term, "Prob_breed")) %>%
  plot_parm_estimates(ms = TRUE)

# Community level effects
draws_beta <- as_draws_df(Occ_mod$beta.comm.samples) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "value")
draws_alpha <- as_draws_df(Occ_mod$alpha.comm.samples) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "value")

plot_half_eye <- function(df){
  df %>% 
    filter(!str_detect(param, "\\.")) %>% # rm .iteration, .draw, and .chain
    ggplot(aes(x = value, y = param, fill = param)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(alpha = .8) +
    labs(x = NULL, y = NULL) +
    guides(fill = "none")
}

# Understand what is contained in (Intercept)
draws_beta %>% 
  #filter(str_detect(param, "Habitat")) %>%
  filter(!str_detect(param, "Ecoregion|Intercept")) %>% 
  plot_half_eye() + 
  labs(title = "100 species")
draws_alpha %>% 
  filter(!str_detect(param, "institucion|Intercept")) %>%
  plot_half_eye() + 
  labs(title = "100 species")

#ggsave("Figures/alpha_post_8fac.png", bg = "white") #beta

# Can play with bayesplot::mcmc_dens() function as well, or as.data.frame() to get posterior draws 

# >z samples --------------------------------------------------------------

# When subsetting the array... rows (posterior samples), columns (species), matrix slice (sites x year_grp)
Occ_mod$z.samples[,1,1] # Species 1 was observed at site 1 

sum(Occ_mod$z.samples[1,,1]) # Species richness at site 1 in posterior sample 1
sum(Occ_mod$z.samples[2,,1]) # Species richness at site 1 in posterior sample 2
sum(Occ_mod$z.samples[1,,2]) # Species richness at site 2 in posterior sample 1
sum(Occ_mod$z.samples[4,,2]) # Species richness at site 2 in posterior sample 4

richness_tbl <- Occ_mod$z.samples %>%
  apply(c(1,3), sum) %>%
  as_tibble()
occ.covs2 <- msOcc_l$occ.covs %>% 
  unite(Id_muestreo_ano, all_of(Row_identifier), sep = ".", remove = FALSE) 
Id_muestreo_ano <- occ.covs2 %>% pull(Id_muestreo_ano)
richness_tbl %>% 
  rename_with(.cols = everything(), ~ Id_muestreo_ano) %>%
  mutate(sample = row_number()) %>%
  pivot_longer(-sample, names_to = "Id_muestreo_ano", values_to = "richness") %>% 
  left_join(occ.covs2[,c("Id_muestreo_ano", "Habitat")]) %>%
  ggplot() +
  geom_density(aes(x = richness, color = Habitat, group = Id_muestreo_ano), #
               alpha = .02)

# >Effective spatial range ------------------------------------------------
# Examine the estimated effective spatial ranges
effective_sp_range <- parms_df %>% 
  filter(parameter == "Theta" & term == "phi") %>% 
  rename(phi_estimate = estimate) %>%
  select(-c(parameter, term)) %>%
  mutate(eff_range_km = (3 / phi_estimate) / 1000)
effective_sp_range %>% pull(eff_range_km) %>% 
  sort() %>%
  setNames(NULL)
