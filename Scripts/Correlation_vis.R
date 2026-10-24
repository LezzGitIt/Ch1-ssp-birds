# Generate correlation plots from predictor variables

# Load libraries and data -------------------------------------------------
library(tidyverse)
library(readxl)
library(xlsx)
library(chron)
library(gridExtra)
library(mvnormtest)
library(conflicted)
library(tidyterra)
library(terra)
library(ggspatial)
ggplot2::theme_set(cowplot::theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)
source("Scripts/spOcc_fns.R")

## Load data
Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"

## NOTE: Removing Hatico points until can digitize
Hatico_pc_nums <- str_pad(1:12, width = 2, pad = 0)

## Landscape metrics (Dist_forest, ssp, te) are built in the sibling Farm_diversity repo by 00a_Extract_lcs.R and 00b_LSM.R
Lsm <- read_lsm()

Event_covs <- read_csv(paste0(Wrangling_repo, Excels, "Event_covs.csv")) %>% 
  filter(!Id_survey_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) %>%
  left_join(Lsm, by = c("Id_survey", "Year"))
Site_covs <- read_csv(paste0(Wrangling_repo, Excels, "Site_covs.csv")) %>% 
  filter(!Id_survey_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 

# >Subset -----------------------------------------------------------------
Forest <- FALSE

if(Forest){
  Site_covs <- Site_covs %>% filter(Habitat == "Bosque")
} else if(!Forest){
  Site_covs <- Site_covs %>% filter(Ecoregion == "Piedemonte")
}

# Join site covs with landscapemetrics. First, match Site_covs_df with 'middle' shapefile, & then overwrite the middle file using the past & ubc files in the correct locations. 
# NOTE:: The 'lc_file' column specifies where the lc information comes from
Site_covs2 <- Event_covs %>% 
  distinct(Id_survey_no_dc, Id_survey, Year, Dist_forest, ssp, te, Canopy_cover, Canopy_height_m) %>% 
  right_join(Site_covs) 

# Calculate correlations --------------------------------------------------
create_cor_matrix <- function(df, cols, cutoff){
  df %>% 
    select({{ cols }}) %>%  
    cor(use = "complete.obs") %>%
    data.frame() %>%
    mutate(across(everything(), round, 2)) %>% 
    mutate(across(everything(), ~ifelse(abs(.x) < cutoff | .x == 1, NA, .x)))
}

# Occupancy covs
Site_covs2 %>%
create_cor_matrix(
  cols = c(Year, Tot_prec, Avg_temp, te, ssp, Canopy_cover, Dist_forest), #Elev
  cutoff = .35
  )

# Detection covariates
Event_covs %>%  
  mutate(Pc_start = as.numeric(Pc_start)) %>%
  select(Julian_day, Pc_start) %>% #Institution_name
  cor(use = "complete.obs") %>%
  data.frame() %>%
  mutate(across(everything(), round, 2))

# Visualize ---------------------------------------------------------------
# Occupancy covs
Site_covs2 %>% mutate(across(where(is.numeric), scale)) %>%
  select(Tot_prec, ssp, Canopy_cover, Year) %>% # Dist_forest
  GGally::ggcorr(label = T, label_size = 2, label_round = 2, hjust = 0.75, size = 3, layout.exp = 1.01)

# Plot correlation between categorical Habitat and continuous predictor variables in causal model 
plot_cat_cor <- function(df, cat_x, y_var, y_title) {
  df %>% distinct({{ cat_x }}, {{ y_var }}) %>%
    ggplot(aes(x = {{ cat_x }}, y = {{ y_var }})) + 
    geom_boxplot(outliers = FALSE) +
    geom_jitter(width = 0.05, alpha = .4) + 
    labs(x = NULL, y = y_title)
}
# Remove Habitat = "Mosaic"? 
Site_covs2 %>% plot_cat_cor(cat_x = Habitat, y_var = ssp, y_title = "Silvopasture (300m)")
Site_covs2 %>% #filter(Habitat != "Mosaic") %>% 
  plot_cat_cor(cat_x = Habitat, y_var = Tot_prec, y_title = "Precipitation (mm)")
Site_covs2 %>% 
  plot_cat_cor(cat_x = Habitat, y_var = Dist_forest, y_title = "Distance to forest (m)")
Site_covs2 %>% plot_cat_cor(cat_x = Habitat, y_var = Canopy_cover, y_title = "Forest in 1km")
Site_covs2 %>% plot_cat_cor(cat_x = Habitat, y_var = Year, y_title = "Year")

# Spatial -----------------------------------------------------------------
## Plot environmental variables to think about spatial autocorrelation
plot_spatial <- function(df, var, Mosaic = FALSE){
  if(!Mosaic){ df  <- df %>% filter(Habitat != "Mosaic") }
    df  %>% terra::vect(geom=c("Long", "Lat")) %>% 
    agitate(maxdist = .013) %>%
    filter(Ecoregion == "Piedemonte") %>%
    ggplot(aes(color = {{ var }})) +
    geom_spatvector(pch = 1, alpha = .8)
}

# Plot
covs_spat_p <- list()
covs_spat_p[[1]] <- plot_spatial(df = Site_covs2, Habitat)
covs_spat_p[[2]] <- plot_spatial(df = Site_covs2, Tot_prec)
covs_spat_p[[3]] <- plot_spatial(df = Site_covs2, Canopy_cover)
covs_spat_p[[4]] <- plot_spatial(df = Site_covs2, Dist_forest)
covs_spat_p[[5]] <- plot_spatial(df = Site_covs2, ssp)
covs_spat_p[[6]] <- plot_spatial(df = Site_covs2, Year)

patchwork::wrap_plots(covs_spat_p, ncols = 3)

# EXTRAS ------------------------------------------------------------------
stop()
Site_covs2 %>% 
  select(where(is.double)) %>%
  GGally::ggpairs()

# Detection covs
Event_covs %>% 
  select(Julian_day, Pc_start, Institution_name) %>%
  GGally::ggpairs()