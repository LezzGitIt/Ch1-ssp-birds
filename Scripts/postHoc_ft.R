## PhD birds in silvopastoral landscapes ##
# Post-hoc functional traits model 

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(posterior)
library(spOccupancy)
library(ggdist)

Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"
Ft <- read_csv(paste0(Wrangling_repo, Excels, "Functional_traits.csv"))

# Occupancy model output
Occ_mod <- readRDS("Rdata/Model_output/Occ_mod_52spp_2026-03-05.rds")

# Forest or non-forest? ---------------------------------------------------
Forest <- TRUE

# Format ------------------------------------------------------------------
# Rows correspond to the posterior MCMC samples and the columns correspond to species
if(Forest){
  variable <- "Canopy_cover"
  y_draws <- as_draws_df(Occ_mod$beta.samples) %>% 
    select(starts_with(variable)) %>% 
    rename_with(~str_remove(., paste0(variable, "-")))
} else{
  variable <- "HabitatSsp"
  y_draws <- as_draws_df(Occ_mod$beta.samples) %>% 
    select(starts_with(variable)) %>% 
    rename_with(~str_remove(., paste0(variable, "-")))
}

# PCA of bill morph
pca <- prcomp(data = Ft, ~ Beak.Depth + Beak.Length_Nares + Beak.Width)
Ft$Beak_size <- pca$x[,1]
Ft$Beak_shape <- pca$x[,2]

# Plot PCA

## General format - Ensure species line up in Ft and y_draws
Spp <- str_replace(names(y_draws), "_", " ")
Ft_spp <- Ft %>%
  filter(Species_ayerbe %in% Spp) %>%
  rename(HWI = `Hand-Wing.Index`) %>%
  mutate(Species_ayerbe = factor(Species_ayerbe, levels = Spp)) %>%
  arrange(Species_ayerbe) %>%
  mutate(across(where(is.numeric), scale)) %>% 
  select(-Species_bl) %>% 
  distinct()

# >Morphology -------------------------------------------------------------
covariates_morph <- c("Mass", "HWI", "Tarsus.Length", "Tail.Length", "Beak_size", "Beak_shape")
Ft_morph <- Ft_spp %>% 
  select(Species_ayerbe, all_of(covariates_morph)) %>%
  # postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(covariates_morph), ~ !is.na(.x)))

# Create list required by postHocLM
postHoc_l_morph <- list(y = y_draws, covs = Ft_morph)

# >Eye --------------------------------------------------------------------
Ft_eye <- Ft_spp %>% filter(!is.na(Eye_resid)) %>% 
  select(-Clutch)
Spp_eye <- Ft_eye %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_eye <- y_draws %>% select(all_of(Spp_eye)) 

# List required by postHocLM
postHoc_l_eye <- list(y = y_draws_eye, covs = Ft_eye)

# >Life-history -----------------------------------------------------------
covariates_lh <- c("Elev_range_final", "Range.Size", "Habitat.Density", "Migration", "Trophic.Level")
Ft_lh <- Ft_spp %>% 
  select(Species_ayerbe, all_of(covariates_lh)) %>%
  # postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(covariates_lh), ~ !is.na(.x)))
postHoc_l_lh <- list(y = y_draws_eye, covs = Ft_lh)

# >Clutch -----------------------------------------------------------------
Ft_clutch <- Ft_spp %>% filter(!is.na(Clutch)) %>% 
  select(-c(Eye_resid, Source_eye))
Spp_clutch <- Ft_clutch %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_clutch <- y_draws %>% select(all_of(Spp_clutch))

# List required by postHocLM
postHoc_l_clutch <- list(y = y_draws_clutch, covs = Ft_clutch)

# Run model ---------------------------------------------------------------
## Morphology models
# Morph variables: Mass + HWI + Tarsus.Length + Tail.Length + Beak_size + Beak_shape
Morph_form <- as.formula(paste("~" , paste(covariates_morph, collapse = " + ")))
postHoc_morph <- postHocLM(formula = Morph_form , data = postHoc_l_morph)
# Residual eye size
postHoc_eye <- postHocLM(formula = ~Eye_resid, data = postHoc_l_eye)

## Life-history models
# Lh vars: Range.Size + Elev_range_final + Habitat.Density + Migration + Trophic.Level
lh_form <- as.formula(paste("~" , paste(covariates_lh, collapse = " + ")))
postHoc_lh <- postHocLM(formula = lh_form, data = postHoc_l_lh)
# Clutch
postHoc_clutch <- postHocLM(formula = ~Clutch, data = postHoc_l_clutch)

# Plot --------------------------------------------------------------------
format.for.plotting <- function(samples){
  samples %>% as_draws_df() %>%
    rename_with(~str_remove(., "Habitat")) %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    mutate(param = as.factor(param))
}

plot_half_eye <- function(df){
  df %>% 
    filter(!param %in% c(".iteration", ".draw", ".chain")) %>%
    ggplot(aes(x = value, y = fct_rev(param), fill = param)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::stat_halfeye(alpha = .8) +
    guides(fill = "none") +
    labs(x = NULL, y = NULL) +
    theme_minimal()
}

## Morphology
format.for.plotting(postHoc_morph$beta.samples) %>% plot_half_eye()
# Residual eye 
format.for.plotting(postHoc_eye$beta.samples) %>% plot_half_eye()

## Life-history
format.for.plotting(postHoc_lh$beta.samples) %>% plot_half_eye()
# Clutch size
format.for.plotting(postHoc_clutch$beta.samples) %>% plot_half_eye()
 