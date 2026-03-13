# Generate correlation plots from predictor variables

# Load libraries and data -------------------------------------------------
library(tidyverse)
library(readxl)
library(xlsx)
library(chron)
library(gridExtra)
library(mvnormtest)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

## Load data
Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"

## NOTE: Removing Hatico points until can digitize
Hatico_pc_nums <- str_pad(1:12, width = 2, pad = 0)

Event_covs <- read_csv(paste0(Wrangling_repo, Excels, "Event_covs.csv")) %>% 
  filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Site_covs <- read_csv(paste0(Wrangling_repo, Excels, "Site_covs.csv")) %>% 
  filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 

# >Subset -----------------------------------------------------------------
Forest <- FALSE

if(!Forest){
  Site_covs <- Site_covs %>% filter(Ecoregion == "Piedemonte")
}

# Join site covs with landscapemetrics. First, match Site_covs_df with 'middle' shapefile, & then overwrite the middle file using the past & ubc files in the correct locations. 
# NOTE:: The 'lc_file' column specifies where the lc information comes from
Site_covs2 <- Event_covs %>% 
  select(Id_muestreo_no_dc, Id_muestreo, Ano, Dist_forest, ssp, te, Canopy_cover, Canopy_height_m) %>% 
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
  cols = c(Ano, Elev, Tot_prec, Avg_temp, te, ssp), 
  cutoff = .35
  )

# Detection covariates
Event_covs %>%  
  mutate(Pc_start = as.numeric(Pc_start)) %>%
  select(Julian_day, Pc_start) %>% #Nombre_institucion
  cor(use = "complete.obs") %>%
  data.frame() %>%
  mutate(across(everything(), round, 2))

# Visualize ---------------------------------------------------------------
# Occupancy covs
GGally::ggcorr(Site_covs2, label = T, label_size = 2, label_round = 2, hjust = 0.75, size = 3, layout.exp = 1.01)

Site_covs2 %>% 
  select(where(is.double)) %>%
  GGally::ggpairs()

# Detection covs
Event_covs %>% 
  select(Julian_day, Pc_start, Nombre_institucion) %>%
  GGally::ggpairs()